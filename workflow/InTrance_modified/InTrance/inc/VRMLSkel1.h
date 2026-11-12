/*!
\file VRMLSkel1.h
\version 1.7
\brief Skeleton VRML code for steering purposes
\author Bernd Dachwald <bernd.dachwald@dlr.de>   (BD)
\author Andreas Ohndorf <andreas.ohndorf@dlr.de> (AO)
\date 18.08.2005 (AO)

\b Purpose:


\b Note:
   This software is protected by national and international copyright.
   Any unauthorized use, reproduction or modificaton is unlawful and
   will be prosecuted. Commercial and non-private application of the
   software in any form is strictly prohibited unless otherwise granted
   by the authors.
   The code is provided without any warranty; without even the implied
   warranty of merchantibility or fitness for a particular purpose.

\b Copyright:
   (c) 2001-2005 Bernd Dachwald, Andreas Ohndorf
*/

#ifndef VRMLSKEL1_H_INCLUDED
#define VRMLSKEL1_H_INCLUDED

// STL
#include <string>

const std::string VRMLSkel1_a("\
#VRML V2.0 utf8\n\
\n\
Background { skyColor 0.0 0.0 0.0 }\n\
\n\
Viewpoint {\n\
  position 0.0 -25.0 0.0\n\
  orientation 1.0 0.0 0.0 1.57\n\
  description \"View 1\"\n\
}\n\
\n\
Viewpoint {\n\
  position 14.434 -14.434 14.434\n\
  orientation 1.0 1.0 0.0 0.85\n\
  description \"View 2\"\n\
}\n\
\n\
Viewpoint {\n\
  position 14.434 +14.434 14.434\n\
  orientation -1.0 1.0 0.0 0.85\n\
  description \"View 3\"\n\
}\n\
\n\
\n\
DEF SC Transform {\n\
 children [\n\
   DEF SC_BODY Transform {\n\
     children [\n\
     Shape {\n\
       appearance Appearance {\n\
         material Material {\n\
           diffuseColor 0.0 0.0 0.95\n\
         }\n\
       }\n\
       geometry Cylinder {\n\
         radius 7.5\n\
         height 0.01\n\
       }\n\
     },\n\
     Shape {\n\
       appearance Appearance {\n\
         material Material {\n\
           diffuseColor 0.5 0.5 0.5\n\
         }\n\
       }\n\
       geometry Box {\n\
         size 0.5 0.5 0.5\n\
       }\n\
     }\n\
     ]\n\
     rotation 0.0 0.0 1.0 -1.57\n\
   },\n\
 \n\
   DEF Arrow_3 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF Arrow_Axis Transform { \n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF Arrow_tip Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                  diffuseColor 0.0 0.0 0.95\n\
                }\n\
             }\n\
             geometry Text {\n");
const std::string VRMLSkel1_b("\
               string \"3\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
     ]\n\
     rotation 1.0 0.0 0.0 1.57\n\
     translation 0.0 0.0 4.5\n\
   },\n\
   DEF Arrow_2 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF Arrow_Axis Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF Arrow_tip Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                  diffuseColor 0.0 0.0 0.95\n\
                }\n\
             }\n\
             geometry Text {\n\
               string \"2\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
     ]\n\
     rotation 1.0 0.0 0.0 0.0\n\
     translation 0.0 4.5 0.0\n\
   },\n\
   DEF Arrow_1 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF Arrow_Axis Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF Arrow_tip Transform {\n\
         children [\n\
           Shape {\n");
const std::string VRMLSkel1_c("\
             appearance Appearance {\n\
               material Material {\n\
                 diffuseColor 0.0 0.0 0.95\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                  diffuseColor 0.0 0.0 0.95\n\
                }\n\
             }\n\
             geometry Text {\n\
               string \"1\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
     ]\n\
     rotation 0.0 0.0 1.0 -1.57\n\
     translation 4.5 0.0 0.0\n\
   }\n\
 ]\n}");


const std::string VRMLSkel1_d("\
# target spacecraft state\n\
#------------------------\n\
\n\
DEF T_SC Transform {\n\
 children [\n\
   DEF T_SC_BODY Transform {\n\
     children [\n\
     Shape {\n\
       appearance Appearance {\n\
         material Material {\n\
           emissiveColor 0.90 0.0 0.0\n\
           transparency 0.5\n\
         }\n\
       }\n\
       geometry Cylinder {\n\
         radius 7.5\n\
         height 0.01\n\
       }\n\
     },\n\
     Shape {\n\
       appearance Appearance {\n\
         material Material {\n\
           emissiveColor 0.5 0.5 0.5\n\
           transparency 0.50\n\
         }\n\
       }\n\
       geometry Box {\n\
         size 0.5 0.5 0.5\n\
       }\n\
     }\n\
     ]\n\
     rotation 0.0 0.0 1.0 -1.57\n\
   },\n\
 \n\
\n\
   DEF T_Arrow_3 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF T_Arrow_Axis Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF T_Arrow_tip Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF T_Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
                }\n\
             }\n\
             geometry Text {\n\
               string \"3\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
\n\
     ]\n\
     rotation 1.0 0.0 0.0 1.57\n");


const std::string VRMLSkel1_e("translation 0.0 0.0 4.5\n\
   },\n\
   DEF T_Arrow_2 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF T_Arrow_Axis Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF T_Arrow_tip Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF T_Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
                }\n\
             }\n\
             geometry Text {\n\
               string \"2\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
\n\
     ]\n\
     rotation 1.0 0.0 0.0 0.0\n\
     translation 0.0 4.5 0.0\n\
   },\n\
   DEF T_Arrow_1 Transform {\n\
     children [\n\
       # arrow body\n\
       DEF T_Arrow_Axis Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cylinder {\n\
               radius 0.1\n\
               height 9.0\n\
             }\n\
            }\n\
          ]\n\
          translation 0.0 0.0 0.0\n\
        },\n\
       # arrow tip\n\
       DEF T_Arrow_tip Transform {\n\
         children [\n");

const std::string VRMLSkel1_f("\
           Shape {\n\
             appearance Appearance {\n\
               material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
               }\n\
             }\n\
             geometry Cone {\n\
               bottomRadius 0.15\n\
               height 1.0\n\
             }\n\
           }\n\
         ]\n\
         translation 0.0 5.0 0.0\n\
       },\n\
       # arrow description\n\
       DEF T_Arrow_desc Transform {\n\
         children [\n\
           Shape {\n\
             appearance Appearance {\n\
                material Material {\n\
                 emissiveColor 0.8 0.1 0.1\n\
                 transparency 0.50\n\
                }\n\
             }\n\
             geometry Text {\n\
               string \"1\"\n\
             }\n\
           }\n\
         ]\n\
         translation 0.5 4.5 0.0\n\
       }\n\
\n\
     ]\n\
     rotation 0.0 0.0 1.0 -1.57\n\
     translation 4.5 0.0 0.0\n\
   }\n\
 ]\n\
}\n");

const std::string VRMLSkel1_g("\
# inertial coordinate system\n\
Shape {\n\
  appearance Appearance{\n\
    material Material { emissiveColor 1.0 1.0 1.0 }\n\
  }\n\
  geometry IndexedLineSet{\n\
    coord Coordinate {\n\
      point [\n\
        -15.0   0.0 0.0,\n\
         15.0   0.0 0.0,\n\
          0.0 -15.0 0.0,\n\
          0.0  15.0 0.0,\n\
          0.0 0.0 -15.0,\n\
          0.0 0.0  15.0\n\
      ]\n\
    }\n\
    coordIndex [\n\
      0, 1, -1,\n\
      2, 3, -1,\n\
      4, 5, -1\n\
    ]\n\
  }\n\
}\n\
\n\
# x axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"X\"\n\
      }\n\
    }\n\
  ]\n\
  translation 15.0 -1.0 -1.0\n\
}\n\
\n\
# y axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"Y\"\n\
      }\n\
    }\n\
  ]\n\
  translation -1.0 15.0 -1.0\n\
}\n\
\n\
# z axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"Z\"\n\
      }\n\
    }\n\
  ]\n\
  translation -1.0 -1.0 15.0\n\
}\n");


 const std::string VRMLSkel1(VRMLSkel1_a+VRMLSkel1_b+VRMLSkel1_c+VRMLSkel1_d
   +VRMLSkel1_e+VRMLSkel1_f+ VRMLSkel1_g);

 

//-----------------------------------------------------------------------------

 // type 2
 
 // viewpoints
 const std::string VRMLSkel2_a("\
#VRML V2.0 utf8\n\n\
Background { skyColor 0.0 0.0 0.0 }\n\n\
Viewpoint {\n\
  position 0.0 -30.0 0.0\n\
  orientation 1.0 0.0 0.0 1.57\n\
  description \"View 1\"\n\
}\n\n\
Viewpoint {\n\
  position 14.434 -14.434 14.434\n\
  orientation 1.0 1.0 0.0 0.85\n\
  description \"View 2\"\n\
}\n\n\
Viewpoint {\n\
  position 14.434 +14.434 14.434\n\
  orientation -1.0 1.0 0.0 0.85\n\
  description \"View 3\"\n\
}\n\n");


const std::string VRMLSkel2_b("\n\
# inertial coordinate system\n\
Shape {\n\
  appearance Appearance{\n\
    material Material { emissiveColor 1.0 1.0 1.0 }\n\
  }\n\
  geometry IndexedLineSet{\n\
    coord Coordinate {\n\
      point [\n\
        -15.0   0.0 0.0,\n\
         15.0   0.0 0.0,\n\
          0.0 -15.0 0.0,\n\
          0.0  15.0 0.0,\n\
          0.0 0.0 -15.0,\n\
          0.0 0.0  15.0\n\
      ]\n\
    }\n\
    coordIndex [\n\
      0, 1, -1,\n\
      2, 3, -1,\n\
      4, 5, -1\n\
    ]\n\
  }\n\
}\n\n");

  // axis labels
const std::string VRMLSkel2_c("\
# x axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"X\"\n\
      }\n\
    }\n\
  ]\n\
  translation 15.0 -1.0 -1.0\n\
}\n\n\
# y axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"Y\"\n\
      }\n\
    }\n\
  ]\n\
  translation -1.0 15.0 -1.0\n\
}\n\n\
# z axis label\n\
Transform {\n\
  children [\n\
    Shape {\n\
      appearance Appearance {\n\
        material Material {\n\
          emissiveColor 1.0 1.0 1.0\n\
        }\n\
      }\n\
      geometry Text {\n\
        string \"Z\"\n\
      }\n\
    }\n\
  ]\n\
  translation -1.0 -1.0 15.0\n\
}\n");


// spacecraft
const std::string VRMLSkel2_d1("\n\
 DEF SC Transform {\n\
 children [\n\
  # body coordinate system\n\
  Transform\n\
  {\n\
   children\n\
   [\n\
    Shape\n\
    {\n\
     appearance Appearance\n\
     {\n\
      material Material { emissiveColor 0.0 0.0 1.0 }\n\
     }\n\
     geometry IndexedLineSet\n\
     {\n\
      coord Coordinate {\n\
       point [\n\
        0.0   0.0 0.0,\n\
        10.0  0.0 0.0,\n\
        0.0   0.0 0.0,\n\
        0.0  10.0 0.0,\n\
        0.0 0.0   0.0,\n\
        0.0 0.0  10.0\n\
       ]\n\
      }\n\
      coordIndex [\n\
       0, 1, -1,\n\
       2, 3, -1,\n\
       4, 5, -1\n\
      ]\n\
     }\n\
    }\n\
   ]\n\
  }\n\
  ,\n");


const std::string VRMLSkel2_d2("\n\
  DEF A1_DSC Transform {\n\
   children [\n\
    Shape {\n\
     appearance Appearance { material Material { diffuseColor 0.0 0.0 0.95 } }\n\
     geometry Text { string \"1\" }\n\
    }\n\
   ]\n\
   rotation 0.0 0.0 1.0 -1.57\n\
   translation 10.0 0.0 0.0\n\
  }\n\
  ,\n\
  DEF A2_DSC Transform {\n\
   children [\n\
    Shape {\n\
     appearance Appearance { material Material { diffuseColor 0.0 0.0 0.95 } }\n\
     geometry Text { string \"2\" }\n\
    }\n\
   ]\n\
   rotation 0.0 0.0 1.0 -1.57\n\
   translation 0.0 10.0 0.0\n\
  }\n  ,\n\
  DEF A3_DSC Transform {\n\
   children [\n\
    Shape {\n\
     appearance Appearance { material Material { diffuseColor 0.0 0.0 0.95 } }\n\
     geometry Text { string \"3\" }\n\
    }\n\
   ]\n\
   rotation 0.0 0.0 1.0 -1.57\n\
   translation 0.0 0.0 10.0\n\
  }\n\
  ,\n\n");

  const std::string VRMLSkel2_d3("\
   DEF SC_BODY Transform {\n\
     children [\n\
     Shape {\n\
      appearance Appearance { material Material { emissiveColor 0.0 0.0 0.90 } }\n\
      geometry Cylinder\n\
      { \n\
        radius 7.5 \n\
        height 0.01\n\
      }\n\
     },\n\
     Shape {\n\
       appearance Appearance { material Material { diffuseColor 0.5 0.5 0.5 } }\n\
       geometry Box {\n\
         size 0.5 0.5 0.5\n\
       }\n\
     }\n\
     ]\n\
     rotation 0.0 0.0 1.0 -1.57\n\
   }\n,\n");

  const std::string VRMLSkel2_d4("\
   # angular velocity\n\
       DEF SC_ANG_VEC Transform\n\
       {\n\
         children \n\
         [\n\
           Transform \n\
           {\n\
             children Shape\n\
             {\n\
               appearance Appearance\n\
               {\n\
                 material Material\n\
                 {\n\
                   emissiveColor 0.0 0.0 0.7 \n\
                   transparency 0.30\n\
                 }\n\
               }\n\
               geometry Cone\n\
               {\n\
                 bottomRadius 0.2\n\
                 height 1.0\n\
               }\n\
             }\n\
             translation 4.5 0.0 0.0\n\
             rotation 0.0 0.0 1.0 -1.57\n\
           },\n\
           Transform \n\
           {\n\
             children Shape\n\
             {\n\
               appearance Appearance { \n\
                   material Material { emissiveColor 0.0 0.0 0.7\n\
                   transparency 0.30\n\
                 }\n\
               }\n\
               geometry Cone\n\
               {\n\
                 bottomRadius 0.2\n\
                 height 1.0\n\
               }\n\
             }\n\
             translation 5.0 0.0 0.0\n\
             rotation 0.0 0.0 1.0 -1.57\n\
           },\n\
           Transform\n\
           {\n\
             children Shape\n\
             {\n\
               appearance Appearance\n\
               {\n\
                 material Material\n\
                 {\n\
                   emissiveColor 0.0 0.0 0.7 \n\
                   transparency 0.30\n\
                 }\n\
               }\n\
               geometry Cylinder\n\
               {\n\
                 radius 0.1\n\
                 height 4.0\n\
               }\n\
             }\n\
             translation 2.0 0.0 0.0\n\
             rotation 0.0 0.0 1.0 -1.57 \n\
           }\n\
         ]\n\
       }\n\
      ]\n}\n");


  const std::string VRMLSkel2_e1("\
# target spacecraft state\n\
#------------------------\n\n\
DEF T_SC Transform {\n\
 children [\n\n\
# body coordinate system\n\
Transform\n\
{ children\n\
  [\n\
    Shape \n\
    {\n\
      appearance Appearance\n\
      {\n\
        material Material { emissiveColor 1.0 0.0 0.0 }\n\
      }\n\
      geometry IndexedLineSet\n\
      {\n\
        coord Coordinate {\n\
          point [\n\
             0.0   0.0 0.0,\n\
            10.0   0.0 0.0,\n\
             0.0   0.0 0.0,\n\
             0.0  10.0 0.0,\n\
             0.0 0.0   0.0,\n\
             0.0 0.0  10.0\n\
         ]\n\
       }\n\
       coordIndex [\n\
         0, 1, -1,\n\
         2, 3, -1,\n\
         4, 5, -1\n\
       ]\n\
     }\n\
   }\n\
  ]\n\
 }\n ,\n");


 const std::string VRMLSkel2_e2("\
 DEF T1_DSC Transform {\n\
  children [\n\
   Shape {\n\
    appearance Appearance { material Material { diffuseColor 0.95 0.0 0.0 } }\n\
    geometry Text { string \"1\" }\n\
   }\n\
  ]\n\
  rotation 0.0 0.0 1.0 -1.57\n\
  translation 10.0 0.0 0.0\n\
 }\n ,\n\
 DEF T2_DSC Transform {\n\
  children [\n\
   Shape {\n\
    appearance Appearance { material Material { diffuseColor 0.95 0.0 0.0 } }\n\
    geometry Text { string \"2\" }\n\
   }\n\
  ]\n\
  rotation 0.0 0.0 1.0 -1.57\n\
  translation 0.0 10.0 0.0\n\
 }\n ,\n\
 DEF T3_DSC Transform {\n\
  children [\n\
   Shape {\n\
    appearance Appearance { material Material { diffuseColor 0.95 0.0 0.0 } }\n\
    geometry Text { string \"3\" }\n\
   }\n\
  ]\n\
  rotation 0.0 0.0 1.0 -1.57\n\
  translation 0.0 0.0 10.0\n\
 }\n ,\n");
 

 const std::string VRMLSkel2_e3("\
 DEF T_SC_BODY Transform {\n\
  children [\n\
   Shape {\n\
    appearance Appearance {\n\
     material Material {\n\
      emissiveColor 0.6 0.0 0.0\n\
      transparency 0.5\n\
     }\n\
    }\n\
    geometry Cylinder {\n\
     radius 7.5\n\
     height 0.01\n\
    }\n\
   },\n\
   Shape {\n\
    appearance Appearance {\n\
     material Material {\n\
      emissiveColor 0.5 0.5 0.5\n\
      transparency 0.50\n\
     }\n\
    }\n\
   geometry Box { size 0.5 0.5 0.5 }\n\
  }\n\
 ]\n\
 rotation 0.0 0.0 1.0 -1.57\n\
}\n ,\n");


const std::string VRMLSkel2_e4("\
 DEF T_SC_ANG_VEC Transform {\n\
  children [\n\
   Transform {\n\
    children Shape {\n\
     appearance Appearance {\n\
      material Material {\n\
       emissiveColor 0.7 0.0 0.0\n\
       transparency 0.30\n\
      }\n\
     }\n\
     geometry Cone {\n\
      bottomRadius 0.2\n\
      height 1.0\n\
     }\n\
    }\n\
    translation 4.5 0.0 0.0\n\
    rotation 0.0 0.0 1.0 -1.57\n\
   },\n\
   Transform {\n\
    children Shape {\n\
     appearance Appearance {\n\
      material Material {\n\
       emissiveColor 0.7 0.0 0.0\n\
       transparency 0.30\n\
      }\n\
    }\n\
      geometry Cone {\n\
      bottomRadius 0.2\n\
      height 1.0\n\
     }\n\
    }\n\
    translation 5.0 0.0 0.0\n\
    rotation 0.0 0.0 1.0 -1.57\n\
   },\n\
   Transform {\n\
    children Shape {\n\
     appearance Appearance {\n\
      material Material {\n\
       emissiveColor 0.7 0.0 0.0\n\
       transparency 0.30\n\
      }\n\
     }\n\
     geometry Cylinder {\n\
       radius 0.1\n\
       height 4.0\n\
     }\n\
    }\n\
    translation 2.0 0.0 0.0\n\
    rotation 0.0 0.0 1.0 -1.57\n\
   }\n\
  ]\n\
  }\n\
 ]\n\
}");


 const std::string VRMLSkel2(VRMLSkel2_a+VRMLSkel2_b+VRMLSkel2_c
   +VRMLSkel2_d1+VRMLSkel2_d2+VRMLSkel2_d3+VRMLSkel2_d4
   +VRMLSkel2_e1+VRMLSkel2_e2+VRMLSkel2_e3+VRMLSkel2_e4);

#endif // VRMLSKEL1_H_INCLUDED