/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    lamvisc

Description
    Calculates the laminar viscosity (mu) from the temperature field T

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <vector>
#include <list>
#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
//    argList::noParallel();
    timeSelector::addOptions(false, true);   // no -constant, with -zeroTime
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();


	volVectorField U
	(
	    IOobject
	    (
		"U",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ
	    ),
	    mesh
	);

        volScalarField T
        (
	     IOobject
             (
    	        "T",
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::MUST_READ
             ), 
    	      mesh
    	);

/*
        surfaceScalarField phi
        (
	     IOobject
             (
    	        "phi",
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::READ_IF_PRESENT
             ), 
    	      linearInterpolate(rho*U) & mesh.Sf()

    	);
*/


        surfaceScalarField phi
        (
	     IOobject
             (
    	        "phi",
    	        runTime.timeName(),
    	        mesh,
    	        IOobject::READ_IF_PRESENT
             ), 
    	      fvc::flux(U)
    	);


    label patchii = mesh.boundaryMesh().findPatchID("tube_to_solid");
    const Foam::fvsPatchField<Foam::Vector<double> > Cpatches = mesh.Cf().boundaryField()[patchii];

    
    scalar nCellsx = Cpatches.size();
    labelListList inBulkRegionList (nCellsx);

//    Info<< "\nCpatches\n" << Cpatches << endl; // Cpatches is the surface centers of the boundary patches.
//    Info<< "\nmesh.C()\n" << mesh.C() << endl; // mesh.c() is the cell centers of the whole mesh.
//    Info<< "\nnCellsx\n" << nCellsx<< endl;    // of coures it's the cell/surface number of boundary patches.
    
    /*
    //return cell centres as volVectorField
    const volVectorField& C() const;

    //return face centres as surfaceVectorField
    const surfaceVectorField& Cf() const;
    
    */
//    const labelList& own = mesh.faceOwner();
  
  
//      labelListList test (nCellsx);
//      labelListList testsum (nCellsx);
//    labelList T_flux_list_for_this_section (nCellsx);
//    labelList flux_list_for_this_section (nCellsx);
//    labelList average_T (nCellsx);
        
        
//      DynamicList<scalarList> test;
//      List<scalarList> test;        
//        DynamicList<float> test;

//        DynamicList<DynamicList<float>> test1;
//        DynamicList<float> test1;        
//        test1.append([1.1]);
//        test1.append([1.1]);
 
        
	DynamicList<DynamicList<float>> T_flux_listlist;  
//	List<List<float>> flux_listlist; 	
                
    forAll(Cpatches, cellj)
    {
//        label nSelectedCell = 1; //if this SelectedCell really necessary? can we just use cellj instead???
                // T_flux_list_for_this_section = [cellj][celli]
                // flux_list_for_this_section = [cellj][celli]
	DynamicList<float> T_flux_list_for_this_section;         
//	List<float> flux_list_for_this_section;

      forAll(mesh.C(), celli)
      {            
          if ( float(Cpatches[cellj].x()) == float(mesh.C()[celli].x()) )  
          {
          
               const labelList& faces = mesh.cells()[celli];

                              
               //find the right face and take it's flux
               // flux is something like CELLi51FACEi1FACE102FLUX4.8805e-09 CELLi58FACEi4FACE17FLUX1.6971e-08
       	forAll(faces, facei)
       	{
       	    if ((mesh.Cf()[faces[facei]].x() - mesh.C()[celli].x()) > 0.0001)
       	    //x difference --> 0.000127397, I hope I picked up a good x difference, could be smaller, check the geometry!
       		{
       		T_flux_list_for_this_section.append(phi[faces[facei]*T[celli]]);
       		Info<<"T_flux_list_for_this_section" << T_flux_list_for_this_section<< endl;  
//       		flux_list_for_this_section.append(phi[faces[facei]]);      		
//       		   label right_face_label = faces[facei];
//       		   Info << "phi[right_face_label]" << phi[right_face_label]<< endl;
//       		   Info << "phi[right_face_label]" << phi[right_face_label]*T[celli]<< endl;
//       		   Info << "T[celli]" << T[celli]<< endl;
//                          test1[cellj].append(phi[right_face_label] *T[celli]);
//                          test.append(phi[right_face_label] *T[celli]);

//       		   Info << "test[cellj]." << test[cellj]<< endl;
//                          test[cellj].setSize(nSelectedCell,phi[right_face_label]*T[celli]);
//       		   Info << "phi[faces[facei]]" << phi[faces[facei]]<< endl;
       		   
//                          test[cellj].setSize(nSelectedCell,phi[faces[facei]]*1e8);
//       		   Info << "test[cellj]" << test[cellj]<< endl;
//       		   T_flux_list_for_this_section[cellj] += (phi[right_face_label] *T[celli]); //just add up wihtout save value
//       		   flux_list_for_this_section[cellj] += (phi[right_face_label]); 

       		}
//       		Info << "CELLi" << celli << "FACEi" << facei << "FACE"<< faces[facei] << "FLUX" << phi[faces[facei]] << "\n"<< endl;
       	}
            
//	       inBulkRegionList[cellj].setSize(nSelectedCell,celli);
//	       nSelectedCell ++;               
          }

      }
      T_flux_listlist.append(T_flux_list_for_this_section);
//      flux_listlist.append(flux_list_for_this_section);      
      
//      testsum[cellj]=sum(test[cellj]);
      
      
      //calculate the average temperature when leaving this loop
//      average_T[cellj] = T_flux_list_for_this_section[cellj] / flux_list_for_this_section[cellj];
    }
//    Info<<"T_flux_listlist" << T_flux_listlist<< endl;   
//    Info<<"flux_listlist" << flux_listlist<< endl;      
    
//    List<float> Avg_T_list;
/*
    List<float> TF;
    List<float> F;    

   
    forAll(T_flux_listlist, i)
//	for (int i=0; i<T_flux_listlist.size(); i++)
    {
    float sum_T_flux = std::accumulate(T_flux_listlist[i].begin(), T_flux_listlist[i].end(), 0.0f);
    float sum_flux = std::accumulate(flux_listlist[i].begin(), flux_listlist[i].end(), 0.0f);
//    Avg_T_list.append(sum_T_flux/sum_flux);
	TF.append(sum_T_flux);
	F.append(sum_flux);

    }
    Info<<"TF" << TF<< endl;   
    Info<<"F" << F<< endl;       
  */
    
    
//     Info<<"test1" << test1<< endl;   
//    Info<<"testsum" << testsum<< endl;
    
//    Info<<"TEST" << test<< endl;
    
//    Info << mesh.cells()<< endl;
// it should return 1)number of faces for the cell, 2)labels for the faces of the cell
// 5(0 1 9850 9950 14950)


//Info << mesh.faceOwner()<< endl;
//should return the cell of certain face
//Info << mesh.faceNeighbour()<< endl;


    
//    Info<< "\nTTTTT\n" << T[0] << endl;
//	    Info<< "\nlabelListList\n" << inBulkRegionList << endl;
//	    Info<< "\nlabelListList\n" << inBulkRegionList[1] << endl;
//	    Info<< "\nlabelListList\n" << inBulkRegionList[10] << endl;	    
/*

        volScalarField T_mean
        (
            IOobject
            (
                "T_mean",
                runTime.timeName(),
                mesh
            ),
        );

        Info<< "Writing: " << "T_mean" << endl;
        T_mean.write();
        
        */
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
