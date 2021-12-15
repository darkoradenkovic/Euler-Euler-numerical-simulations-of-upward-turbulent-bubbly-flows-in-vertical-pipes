/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

scale 1;
   
define(cos45, 0.70710678)
define(cos60, 0.5)
define(sin60, 0.866025404)

define(cos30, 0.866025404)
define(sin30, 0.5)

define(sin225, 0.382683432)
define(cos225, 0.923879533)

define(R, calc(0.025/2))
define(Rx, calc(R*cos45))   // rux=ruz
define(r1, calc(0.7*R))
define(r1x, calc(r1*cos45))
define(r2, calc(0.8*r1))
define(r2x, calc(r2*cos45))

define(X, calc(r2*0.973214))

define(LA, 0.0)
define(LB, 2)


// // // Brojevi podeoka u blokovima

define(b, 20) // za unutrasnji blok 0, ovaj broj oznacava i broj podeoka radijalnim zracima!!!

define(kp, 20) // broj kruznih podeoka u blokovima 1 2 
define(sgKP, 0.15) // simpleGrading za blokove 1 2 3 4 10 11 12 13 ... u z pravcu (u lokalnim
                  // koordinatnim sistemima za svaki blok z osa je postavljena u radijalnom pravcu)

define(LABy, 100)  


vertices
(
 // RavanA
  (0.0 LA 0.0) vlabel(A0) 
  (r2 LA 0) vlabel(A1)
  (R LA 0) vlabel(A2)
  (Rx LA -Rx) vlabel(A3)
  (r1x LA -r1x) vlabel(A4) 
  (0 LA -r2) vlabel(A5)
  (0 LA -R) vlabel(A6)

 // RavanB
  (0.0 LB 0.0) vlabel(B0) 
  (r2 LB 0) vlabel(B1)
  (R LB 0) vlabel(B2)
  (Rx LB -Rx) vlabel(B3)
  (r1x LB -r1x) vlabel(B4) 
  (0 LB -r2) vlabel(B5)
  (0 LB -R) vlabel(B6)


);

blocks
(    
	// Svaki blok ima svoj koordinatni sistem !

   hex (A0 A1 A4 A5 B0 B1 B4 B5) (b b LABy) 
   simpleGrading
	(
	   1
	   1
	   1
	)						   //   blok 0

   hex (A1 A2 A3 A4 B1 B2 B3 B4) (kp b LABy) 
   simpleGrading
	(
	   sgKP
	   1
	   1
	)						   //   blok 1

   hex (A4 A3 A6 A5 B4 B3 B6 B5) (kp b LABy) 
   simpleGrading
	(
	   sgKP
	   1
	   1
	)						   //   blok 2

 


);

  edges
  (
   	// Ravan A
	arc A2 A3 (calc(R*cos225) LA calc(-R*sin225))
	arc B2 B3 (calc(R*cos225) LB calc(-R*sin225))

	arc A3 A6 (calc(R*cos60) LA calc(-R*sin60))
	arc B3 B6 (calc(R*cos60) LB calc(-R*sin60))

	arc A1 A4 (X LA calc(-r1x/2.0))
	arc B1 B4 (X LB calc(-r1x/2.0))

	arc A5 A4 (calc(r1x/2.0) LA -X)
	arc B5 B4 (calc(r1x/2.0) LB -X)

  );

boundary
(
    inlet
    {
       type patch;
       faces
       (
         (A0 A5 A4 A1)
         (A1 A4 A3 A2)
         (A6 A3 A4 A5)
      );
    }

    outlet
    {
      type patch;
      faces
      (
        (B0 B1 B4 B5)
        (B1 B2 B3 B4)
        (B5 B4 B3 B6)
      );
    } 

	walls
	{
	  type wall;
	  faces
	  (
	    (A2 A3 B3 B2)
	    (A3 A6 B6 B3)
	  );
	} 

	symmetry1
	{
	  type symmetryPlane;
	  faces
	  (
	    (A0 A1 B1 B0)
	    (A1 A2 B2 B1)
	  );
	} 
	symmetry2
	{
	  type symmetryPlane;
	  faces
	  (
	    (A0 B0 B5 A5)
	    (A6 A5 B5 B6)
	  );
	} 
);

mergePatchPairs
(
);
    
    







