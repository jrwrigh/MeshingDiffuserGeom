#!/usr/bin/env python
import math
from pathlib import Path
import argparse
import OCC
import OCC.gp as gp
import OCC.GC as GC
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.BRepPrimAPI import BRepPrimAPI_MakeRevol
from OCC.Display.SimpleGui import init_display
from OCC.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.StlAPI import StlAPI_Writer
from OCC.BRepMesh import BRepMesh_IncrementalMesh

from OCC.SMESH import SMESH_Gen, SMESH_MeshVSLink
from OCC.StdMeshers import (StdMeshers_Arithmetic1D,
                            StdMeshers_TrianglePreference,
                            StdMeshers_Regular_1D, StdMeshers_Quadrangle_2D,
                            StdMeshers_MEFISTO_2D, StdMeshers_MaxLength)

from OCC.GProp import GProp_GProps
from OCC.BRepGProp import brepgprop_SurfaceProperties

import OCCUtils

from OCC.TopoDS import TopoDS_Face

# Here is an example script of using pythonocc to create a conical diffuser
# from nothing and save it as a STEP file. The diffuser is created via a
# revolve of the initial geometry.

# ---Visualization of the diffuser geometry Note that the shoulder of the
# diffuser (marked by "p2, p3, p4, e6") has those three points in that general
# order. Point2 has the exact same height as pnt1 while pnt3 and pnt4 form the
# rest of the fillet arc. That fillet arc is edge 6.
#                                                   e3
#                                     p5*---------------------------*p6
#                                      /                            |
#                        e1           /e2                           |
#               p1*-----------------*/                              |
#                 |                 ^p2,p3,p4,e6                    |e4
#                 |                                                 |
#               e0|                                                 |
#                 |                                                 |
#                 |                                                 |
# rot_axis-->   p0*-------------------------------------------------*p7
#                                            e5
# Used the API Documentation for the pythonocc module quite a bit as a general reference:
#   #   https://cdn.rawgit.com/tpaviot/pythonocc-core/804f7f3/doc/apidoc/0.18.1/OCC.STEPControl.html?highlight=stepcontrol
# Mostly followed the pythonocc script example bottle.py found here:
#    #    https://gist.github.com/jf---/7f66d8e711d7a1cac66684178b2e7740
# Also used Open CASCADE Tutorial which is a more guided version of the bottle.py script, but doesn't use the pythonocc library.
#   # https://www.opencascade.com/doc/occt-7.3.0/overview/html/occt__tutorial.html
# Also used for the STEP writing procedure was this OpenCASCADE resource:
#   # https://www.opencascade.com/doc/occt-7.3.0/overview/html/occt_user_guides__step.html#occt_step_3_1
# Note: pythonocc is based on the Open CASCADE stuff

#############################################
#-----------argparse

# Test command line input:
# python ConicalDiffuserGeom.py -l 200 300 -r 20 30 -t 20 -a 10 -f test.stl --filetype stl -v

parser = argparse.ArgumentParser(
    description='Create conical diffuser geometry in a STEP file.')

parser.add_argument(
    '-l',
    '--lengths',
    dest='l',
    metavar=('IN_L', 'OUT_L'),
    type=float,
    nargs=2,
    help='The lengths of the the inlet and outlet, respectively')
parser.add_argument(
    '-r',
    '--radius',
    dest='r',
    metavar=('IN_R', 'OUT_R'),
    type=float,
    nargs=2,
    help='The radius of the inlet and outlet, respectively')
parser.add_argument(
    '-a',
    '--angle',
    dest='a',
    type=float,
    help='Half Angle of the diffuser in degrees, phi')
parser.add_argument(
    '-t',
    '--transitionradius',
    dest='t',
    type=float,
    help=
    'The radius of the transition from the inlet radius to the diffuser angle')
parser.add_argument(
    '-f',
    '--file',
    dest='f',
    type=str,
    help='The path or file name that the geometry will be saved to.')
parser.add_argument('--TEST', action="store_true")
parser.add_argument('-v', '--verbose', dest='v', action='count', default=0)
parser.add_argument(
    '--filetype',
    type=str,
    help='The file type that should be output (stl, STEP, etc.)')

args = parser.parse_args()

if not args.TEST:
    inletr = args.r[0]
    inletl = args.l[0]
    phi = args.a  # in degrees
    outletr = args.r[1]
    outletl = args.l[1]
    filletr = args.t
    filepath = Path(args.f).resolve()
    filetype = args.filetype

elif args.TEST:
    inletr = 50
    inletl = 200
    phi = 10  # in degrees
    outletr = 60
    outletl = 300
    filletr = 20
    filepath = Path('testingSMESH.stl').resolve()
    filetype = 'stl'
    makeDiffuser_kwargs = {
        'radii': (inletr, outletr),
        'lengths': (inletl, outletl),
        'phi': phi,
        'transitionr': filletr
    }

if args.v >= 1:
    inputs = {
        'Inlet Radius:': inletr,
        'Inlet Length': inletl,
        'Half Angle:': phi,
        'Outlet Radius:': outletr,
        'Outlet Length:': outletl,
        'Transition Radius': filletr,
        'File Path:': filepath,
        'File Type:': filetype
    }

    print('INPUTS:\n-------')
    for key in inputs.keys():
        if type(inputs[key]) == int or type(inputs[key]) == float:
            print('{:16}\t{:10.5f}'.format(key, inputs[key]))
        elif type(inputs[key]) == str:
            print('{:16}\t{:10}'.format(key, inputs[key]))

    calculated_values = {
        'Area Ratio:': outletr**2 / inletr**2,
        'Diffuser Length:': (outletr - inletr) / math.tan(math.radians(phi))
    }
    print('\nDIFFUSER VALUES:')
    print('----------------')
    for key in calculated_values.keys():
        if type(calculated_values[key]) == int or type(
                calculated_values[key]) == float:
            print('{:16}\t{:10.5f}'.format(key, calculated_values[key]))
        elif type(calculated_values[key]) == str:
            print('{:16}\t{:10}'.format(key, calculated_values[key]))

    print('\n')

#############################################
#---------Geometry Creation Functions
#############################################


def makePoints(radii, lengths, transitionr, phi):
    inletr, outletr = radii
    inletl, outletl = lengths
    transitionr = filletr
    # Initial Parameter calculations
    phi = math.radians(phi)
    pntspy = []
    # This has the math that defines the three points (p2, p3, p4) that lie on an arc.
    pntspy.append((0, 0, 0))
    pntspy.append((0, inletr, 0))
    pntspy.append((inletl, inletr, 0))
    pntspy.append((filletr * math.sin(phi / 2) + pntspy[2][0],
                   filletr * (1 - math.cos(phi / 2)) + pntspy[2][1], 0))
    pntspy.append((filletr * math.sin(phi) + pntspy[2][0],
                   filletr * (1 - math.cos(phi)) + pntspy[2][1], 0))
    pntspy.append(((outletr - inletr - (pntspy[4][1] - inletr)) / math.tan(phi)
                   + pntspy[4][0], outletr, 0))
    pntspy.append((pntspy[5][0] + outletl, outletr, 0))
    pntspy.append((pntspy[5][0] + outletl, 0, 0))

    if args.v >= 1: print('Points are calculated')
    # Creating the base sketch points
    pnts = []
    for pnt in pntspy:
        pnts.append(OCC.gp.gp_Pnt(*pnt))
    if args.v >= 1: print('Points are created')

    return pnts


def makeEdges(pnts):
    """ Creating the base curves and segments"""

    # Segmentkeys shows what points form an edge (as opposed to the fillet arc)
    segmentkeys = [(0, 1), (1, 2), (4, 5), (5, 6), (6, 7), (7, 0)]
    segments = []
    for key in segmentkeys:
        segments.append(OCC.GC.GC_MakeSegment(pnts[key[0]], pnts[key[1]]))
    # Adds edge6, which is the fillet arc
    segments.append(OCC.GC.GC_MakeArcOfCircle(pnts[2], pnts[3], pnts[4]))
    # Creating edges. Could also do this directly from points
    edges = []
    for seg in segments:
        edges.append(BRepBuilderAPI_MakeEdge(seg.Value()))
    if args.v >= 1: print('Edges are created')
    return edges


def makeWire(edges):
    """ Creating Wire
    
    Note the edges have to be added in an order than default. I chose the keep
    the edge numbering in a sensible format (the format shown in the diffuser
    visualization). However, wire.Add() requires that each edge be connected to
    the current wire. The current order goes around the physical loop starting at
    edge0
    """

    wire = BRepBuilderAPI_MakeWire()
    edgeorder = (0, 1, 6, 2, 3, 4, 5)
    for edgen in edgeorder:
        wire.Add(edges[edgen].Edge())
    if args.v >= 1: print('Wire are created')
    return wire


def makeFacefromWire(wire):
    """ Create face from Wire"""

    face = BRepBuilderAPI_MakeFace(wire.Wire())
    if args.v >= 1: print('Face are created')

    return face


def makeRevolve(face):
    # Get X axis, which face will be revolved around
    xaxis = OCC.gp.gp_OX()

    # Revolve the face around the xaxis
    diffuser = BRepPrimAPI_MakeRevol(face.Face(), xaxis)
    if args.v >= 1: print('Geometry is created')
    return (diffuser)


def makeDiffuser(radii, lengths, transitionr, phi):
    pnts = makePoints(radii, lengths, transitionr, phi)
    edges = makeEdges(pnts)
    wire = makeWire(edges)
    face = makeFacefromWire(wire)
    diffuser = makeRevolve(face)
    return diffuser


###################################
#-----------Meshing Functions
###################################


def makeIncrementalMesh(shape, lindefl=0.1, angdefl=.07):
    """ Make mesh using BRepMesh_IncrementalMesh 
    
     BRepMesh_IncrementalMesh was found to have some issues with the
     transitions radius (basically turning it into a chamfer instead,
     getting rid of the actual curvature), so SMESH was then explored.

     lindefl = Linear Deflection Setting
     angdefl = Angular Deflection Setting
    """

    try:
        mesh = BRepMesh_IncrementalMesh(diffuser.Shape(), lindefl, True,
                                        angdefl, False)
    except:
        print('Mesh failed')
    else:
        if args.v >= 1:
            print(f'\nDiffuser shape has been meshed.')
    return mesh


def makeSMESH(shape):
    """ Make mesh using SMESH

    Most of this was taken from core_mesh_surfacic.py in pythonocc-demos
    """
    aMeshGen = SMESH_Gen()
    aMesh = aMeshGen.CreateMesh(0, True)
    
    # 1D
    an1DHypothesis = StdMeshers_MaxLength(0, 0, aMeshGen)
    an1DHypothesis.SetLength(10)
    an1DAlgo = StdMeshers_Regular_1D(1, 0, aMeshGen)  # interpolation
    # 2D
    a2dHypothseis = StdMeshers_TrianglePreference(
        2, 0, aMeshGen)  # define the boundary
    a2dAlgo = StdMeshers_MEFISTO_2D(3, 0, aMeshGen)
    # alculate mesh
    aMesh.ShapeToMesh(diffuser.Shape())
    # Assign hyptothesis to mesh
    aMesh.AddHypothesis(diffuser.Shape(), 0)
    aMesh.AddHypothesis(diffuser.Shape(), 1)
    aMesh.AddHypothesis(diffuser.Shape(), 2)
    aMesh.AddHypothesis(diffuser.Shape(), 3)
    # Compute the data
    aMeshGen.Compute(aMesh, aMesh.GetShapeToMesh())
    return aMesh


#######################################################
#----------------Writing the Geometry to CAD
#######################################################


def writeSTEP(filepath, shape):
    """ Write STEP file of the shape"""

    # STEPControl_AsIs says to make the STEP model the same geometry type as the
    # shape (ie. a solid Shape should be a STEP Solid)
    writer = STEPControl_Writer()
    writer.Transfer(diffuser.Shape(), STEPControl_AsIs)

    # Write the STEP file to the given filepath
    try:
        writer.Write(filepath.as_posix())
    except:
        print('Write failed')
    else:
        if args.v >= 1:
            print(f'\nSTEP file was successfully saved to:\n {filepath}')


def writestl(filepath, mesh):
    """ Write stl of IncrementalMesh 

    Note: this is for IncrementalMesh files only. SMESH files have their own
    built in `.ExportStl()` method for this purpose 
    """

    writer = StlAPI_Writer()
    writer.SetASCIIMode(True)

    try:
        writer.Write(mesh.Shape(), filepath.as_posix())
    except:
        print('Write failed')
    else:
        if args.v >= 1:
            print(f'\nstl file was successfully saved to:\n {filepath}')


if args.TEST:
    diffuser = makeDiffuser(**makeDiffuser_kwargs)

####################################
#-----------Misc Geo Functions
####################################


def getFaceArea(face):
    """ Return the area of the face object given """

    system = GProp_GProps()
    brepgprop_SurfaceProperties(face, system)
    return (system.Mass())


def returnSmallestFace(shape):
    """ Return ID# of the smallest face in the shape"""

    topology = OCCUtils.Topo(shape)
    for n, face in enumerate(topology.faces()):
        assert isinstance(face, TopoDS_Face)
        facearea = getFaceArea(face)
        if n == 0:
            maxarea = facearea
            faceobj = face
        elif facearea < maxarea:
            maxarea = facearea
            faceobj = face

    return face


#####################################
#--------------Visualization
#####################################


def visualize(Shapeobject=diffuser):
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayColoredShape(Shapeobject.Shape())

    start_display()
