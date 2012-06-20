#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV, Andrey SIMURZIN
##


#---------------------------------------------------------------------------
from Foam import ref, man


#---------------------------------------------------------------------------
def _createFields( runTime, mesh, simple ):
    ref.ext_Info() << "Reading thermophysical properties\n" << ref.nl
    
    thermo = man.basicPsiThermo.New( mesh )

    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.READ_IF_PRESENT,
                                            ref.IOobject.AUTO_WRITE ),
                              man.volScalarField( thermo.rho(), man.Deps( thermo ) ) )

    p = man.volScalarField( thermo.p(), man.Deps( thermo ) )
    h = man.volScalarField( thermo.h(), man.Deps( thermo ) )
    psi = man.volScalarField( thermo.psi(), man.Deps( thermo ) )
   
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    phi = man.compressibleCreatePhi( runTime, mesh, rho, U )
    
    pRefCell = 0
    pRefValue = 0.0
    
    pRefCell, pRefValue = ref.setRefCell( p, simple.dict(), pRefCell, pRefValue )
    
    rhoMax = ref.dimensionedScalar( simple.dict().lookup( ref.word( "rhoMax" ) ) )
    rhoMin = ref.dimensionedScalar( simple.dict().lookup( ref.word( "rhoMin" ) ) )
    
    ref.ext_Info() << "Creating turbulence model\n" << ref.nl
    turbulence = man.compressible.RASModel.New( rho,
                                                U,
                                                phi,
                                                thermo )
    
    initialMass = ref.fvc.domainIntegrate( rho )
    
    return thermo, rho, p, h, psi, U, phi, pRefCell, pRefValue, turbulence, initialMass, rhoMax, rhoMin


#--------------------------------------------------------------------------------------
def _UEqn( phi, U, p, rho, turbulence ):
    # Solve the Momentum equation
    UEqn = man.fvVectorMatrix( turbulence.divDevRhoReff( U ), man.Deps( turbulence ) )  + man.fvm.div( phi, U )
    
    UEqn.relax()
    
    ref.solve( UEqn == -ref.fvc.grad( p ) )
    
    return UEqn


#--------------------------------------------------------------------------------------
def _hEqn( U, phi, h, turbulence, rho, p, thermo ):
    
    hEqn = ref.fvm.div( phi, h ) - ref.fvm.Sp( ref.fvc.div( phi ), h ) - ref.fvm.laplacian( turbulence.alphaEff(), h ) \
           == - ref.fvc.div( phi, 0.5 * U.magSqr(), ref.word( "div(phi,K)" ) )

    hEqn.relax()
    hEqn.solve()

    thermo.correct()
    pass


#--------------------------------------------------------------------------------------
def _pEqn( runTime,mesh, UEqn, rho, thermo, psi, U, p, phi, \
           pRefCell, pRefValue, cumulativeContErr, simple, rhoMin, rhoMax  ):
   
    rho << thermo.rho()
    rho << rho.ext_max( rhoMin )
    rho << rho.ext_min( rhoMax )
    rho.relax()

    p0 = ref.volScalarField(p)
    
    AU = ref.volScalarField( UEqn.A() )
    AtU = AU - UEqn.H1()
    U << UEqn.H() / AU
    
    # UEqn.clear()
    
    closedVolume = False
    
    if simple.transonic(): 
       while simple.correctNonOrthogonal():
           phid = ref.surfaceScalarField( ref.word( "phid" ),
                                          ref.fvc.interpolate( psi * U ) & mesh.Sf() )

           phic = ref.surfaceScalarField( ref.word( "phic" ),
                                          ref.fvc.interpolate( rho() / AtU - rho() / AU ) * ref.fvc.snGrad( p ) * mesh.magSf() + \
                                          phid * ( ref.fvc.interpolate( p ) - ref.fvc.interpolate( p, ref.word( "UD" ) ) ) ) # mixed calcutions

           #refCast<mixedFvPatchScalarField>(p.boundaryField()[1]).refValue()
           #    = p.boundaryField()[1];
           
           pEqn = ref.fvm.div( phid, p ) + ref.fvc.div( phic ) - ref.fvm.Sp( ref.fvc.div( phid ), p ) + \
                  ref.fvc.div( phid ) * p - ref.fvm.laplacian( rho() / AtU, p ) # mixed calcutions
           # pEqn.relax();

           pEqn.setReference( pRefCell, pRefValue )
           pEqn.solve()

           if simple.finalNonOrthogonalIter():
               phi == phic + pEqn.flux()
               pass
           pass
    else:
        
        while simple.correctNonOrthogonal():
            phi << ( ref.fvc.interpolate( rho * U ) & mesh.Sf() )
            closedVolume = ref.adjustPhi( phi, U, p )
            
            phi += ref.fvc.interpolate( rho / AtU - rho / AU ) * ref.fvc.snGrad( p ) * mesh.magSf()
            pEqn = ref.fvc.div( phi ) - ref.fvm.laplacian( rho / AtU, p )
            pEqn.setReference( pRefCell, pRefValue )
            
            pEqn.solve()
            if simple.finalNonOrthogonalIter():
                phi += pEqn.flux()
                pass
            pass
            
    # The incompressibe for of the continuity error check is appropriate for
    # steady-state compressible also.
    cumulativeContErr = ref.ContinuityErrs( phi, runTime, mesh, cumulativeContErr )
    
    # Explicitly relax pressure for momentum corrector
    p.relax()
    
    U -= ( ref.fvc.grad( p0 ) * ( 1.0 / AU - 1.0 / AtU ) + ref.fvc.grad( p ) / AtU )
    # U -= fvc::grad(p)/AU;
    
    U.correctBoundaryConditions()
    
    # For closed-volume cases adjust the pressure and density levels
    # to obey overall mass continuity
    if closedVolume:
        p += ( initialMass - ref.fvc.domainIntegrate( psi * p ) ) / ref.fvc.domainIntegrate( psi )
        pass
    
    rho << thermo.rho()
    rho << rho.ext_max( rhoMin )
    rho << rho.ext_min( rhoMax )
    
    if not simple.transonic():
        rho.relax()
        pass

    ref.ext_Info() << "rho max/min : " << rho.ext_max().value() << " " << rho.ext_min().value() << ref.nl
    
    pass
    
    return cumulativeContErr


#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    simple = man.simpleControl( mesh )

    thermo, rho, p, h, psi, U, phi, pRefCell, pRefValue, turbulence, initialMass, rhoMax, rhoMin = _createFields( runTime, mesh, simple )
    
    cumulativeContErr = ref.initContinuityErrs()

    ref.ext_Info() << "\nStarting time loop\n" << ref.nl
    
    while simple.loop() :
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl
                
        # Pressure-velocity SIMPLE corrector
        UEqn = _UEqn( phi, U, p, rho, turbulence )
            
        cumulativeContErr = _pEqn( runTime,mesh, UEqn, rho, thermo, psi, U, p, phi, \
                                                             pRefCell, pRefValue, cumulativeContErr, simple, rhoMin, rhoMax  )

        _hEqn( U, phi, h, turbulence, rho, p, thermo )
        
        
        turbulence.correct()
        
        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n"

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020101" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass   
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.1 or higher \n "


    
#--------------------------------------------------------------------------------------
