""" Minimal PDK for EBeam constructed with ZeroPDK. """

import os
import sys
import logging
from collections import abc
from zeropdk import Tech
from zeropdk.pcell import PCell
import klayout.db as kdb

logger = logging.getLogger()
lyp_path = os.path.join(os.path.dirname(__file__), "EBeam.lyp")

pwd = os.path.dirname(os.path.realpath(__file__))
sys.path.append(pwd)

# Technology file
EBeam = Tech.load_from_xml(lyp_path)
# Currently, EBeam contains only one attribute:
# EBeam.layers, to be used in the pcell definitions
# below

TECHLAYERS = EBeam.layers

# PCells
from zeropdk.layout.geometry import rotate, rotate90
from zeropdk.pcell import (
    TypeDouble,
    TypeInt,
    TypeBoolean,
    TypeString,
    TypeList,
    TypePoint,
    TypeVector,
    TypeLayer,
)

from zeropdk.pcell import PCellParameter, ParamContainer


def define_param(name, type, description, default=None, **kwargs):
    from zeropdk.pcell import PCellParameter

    return PCellParameter(
        name=name, type=type, description=description, default=default, **kwargs
    )


from zeropdk.pcell import objectview


class EBeamLayersMixin(PCell):
    """ Abstract class with more concise layer handling

    """

    params = ParamContainer(
        define_param("silayer", TypeLayer, "Si (1/0)", default=TECHLAYERS["Si"]),
        define_param(
            "31_Si_p6nm",
            TypeLayer,
            "'31_Si_p6nm' (31/0)",
            default=TECHLAYERS["31_Si_p6nm"],
        ),
        define_param("textl", TypeLayer, "Text (10/0)", default=TECHLAYERS["Text"]),
        define_param("si_n", TypeLayer, "'Si N' (20/0)", default=TECHLAYERS["Si N"]),
        define_param(
            "si_npp", TypeLayer, "'Si N++' (24/0)", default=TECHLAYERS["Si N++"]
        ),
        define_param("SEM", TypeLayer, "SEM (200/0)", default=TECHLAYERS["SEM"]),
        define_param("M1", TypeLayer, "M1 (41/0)", default=TECHLAYERS["M1"]),
        define_param("ML", TypeLayer, "'12_M2' (12/0)", default=TECHLAYERS["12_M2"]),
        define_param(
            "ML_Open", TypeLayer, "'13_MLopen' (13/0)", default=TECHLAYERS["13_MLopen"]
        ),
        define_param("VC", TypeLayer, "VC (40/0)", default=TECHLAYERS["VC"]),
        define_param(
            "M_Heater", TypeLayer, "'M Heater' (47/0)", default=TECHLAYERS["M Heater"]
        ),
        define_param(
            "devrec", TypeLayer, "DevRec (68/0)", default=TECHLAYERS["DevRec"]
        ),
        define_param(
            "pinrec", TypeLayer, "PinRec (1/10)", default=TECHLAYERS["PinRec"]
        ),
        define_param(
            "FbrTgt", TypeLayer, "FbrTgt (81/0)", default=TECHLAYERS["FbrTgt"]
        ),
        define_param(
            "Lumerical", TypeLayer, "Lumerical (733/0)", default=TECHLAYERS["Lumerical"]
        ),
        define_param(
            "MEEP_SOURCE1", TypeLayer, "MEEP_SOURCE1 (1001/0)"
        ),
        define_param(
            "MEEP_SOURCE2", TypeLayer, "MEEP_SOURCE1 (1002/0)"
        ),
        define_param(
            "MEEP_PORT1", TypeLayer, "MEEP_PORT1 (1011/0)"
        ),
        define_param(
            "MEEP_PORT2", TypeLayer, "MEEP_PORT2 (1012/0)"
        ),
        define_param(
            "MEEP_PORT3", TypeLayer, "MEEP_PORT3 (1013/0)"
        ),
        define_param(
            "MEEP_PORT4", TypeLayer, "MEEP_PORT4 (1014/0)"
        )
    )

    def get_layers(self, layout):
        """ Instantiates layers into layout and returns a convenient objectview.

            Use with `lay = self.get_layers(layout); lay.Si`

        """
        cp = self.get_cell_params()
        # breakpoint()

        lay = objectview({})
        lay.Si = layout.layer(cp.silayer)
        lay.SiN = layout.layer(cp.si_n)
        lay.SiNpp = layout.layer(cp.si_npp)
        lay.M1 = layout.layer(cp.M1)
        lay.ML = layout.layer(cp.ML)
        lay.M_Heater = layout.layer(cp.M_Heater)
        lay.PinRec = layout.layer(cp.pinrec)
        lay.DevRec = layout.layer(cp.devrec)
        lay.Text = layout.layer(cp.textl)
        lay.MLOpen = layout.layer(cp.ML_Open)
        
        lay.MEEP_SOURCE1 = layout.layer(cp.MEEP_SOURCE1)
        lay.MEEP_SOURCE2 = layout.layer(cp.MEEP_SOURCE2)
        lay.MEEP_PORT1 = layout.layer(cp.MEEP_PORT1)
        lay.MEEP_PORT2 = layout.layer(cp.MEEP_PORT2)
        lay.MEEP_PORT3 = layout.layer(cp.MEEP_PORT3)
        lay.MEEP_PORT4 = layout.layer(cp.MEEP_PORT4)

        return lay


class PositionMixin(PCell):
    """ handles the angle_ex parameter """

    params = ParamContainer(
        PCellParameter(
            name="angle_ex",
            type=TypeDouble,
            description="Placement Angle (0, 90, ..)",
            default=0,
        )
    )

    def origin_ex_ey(self, multiple_of_90=False):  # pylint: disable=unused-argument
        EX = kdb.DVector(1, 0)
        cp = self.get_cell_params()
        origin = kdb.DPoint(0, 0)
        # if 'angle_ex' not in cp.__dict__:
        #     cp.angle_ex = 0
        if multiple_of_90:
            if cp.angle_ex % 90 != 0:
                raise RuntimeError("Specify an angle multiple of 90 degrees")

        from math import pi

        ex = rotate(EX, cp.angle_ex * pi / 180)
        ey = rotate90(ex)
        return origin, ex, ey


class EBeamCell(EBeamLayersMixin, PositionMixin, PCell):
    """ Main layout class for instantiating PDK-compatible cells """

    pass


PDKCell = EBeamCell

# set LOCAL_GDS_DIR as (location of ebeam_pdk.py)/gds_cells
LOCAL_GDS_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gds_cells")


def GDSCell(cell_name, filename=None, gds_dir=LOCAL_GDS_DIR):
    """
        Args:
            cell_name: cell within that file.
            filename: is the gds file name.
        Returns:
            (class) a GDS_cell_base class that can be inherited
    """
    from zeropdk.pcell import GDSCell as DefaultGDSCell

    defaultGDSClass = DefaultGDSCell(cell_name, filename=filename, gds_dir=gds_dir)

    class pdkGDSClass(defaultGDSClass, PDKCell):
        def draw_gds_cell(self, cell):
            layout = cell.layout()
            gdscell = self.get_gds_cell(layout)

            origin, _, _ = self.origin_ex_ey()
            cell.insert_cell(gdscell, origin, self.params.angle_ex)
            return cell

        def draw(self, cell):
            return self.draw_gds_cell(cell), {}

    return pdkGDSClass


# Overriding default layers

from zeropdk.default_library.io import DCPad, DCPadArray


class DCPad(DCPad):
    params = ParamContainer(
        PCellParameter(
            name="layer_metal",
            type=TypeLayer,
            description="Metal Layer",
            default=EBeam.layers["M1"],
        ),
        PCellParameter(
            name="layer_opening",
            type=TypeLayer,
            description="Open Layer",
            default=EBeam.layers["13_MLopen"],
        ),
    )


class DCPadArray(DCPadArray):
    params = ParamContainer(
        PCellParameter(
            name="layer_metal",
            type=TypeLayer,
            description="Metal Layer",
            default=EBeam.layers["M1"],
        ),
        PCellParameter(
            name="layer_opening",
            type=TypeLayer,
            description="Open Layer",
            default=EBeam.layers["13_MLopen"],
        ),
    )


# Helper functions

WAVEGUIDE_RADIUS = 10
WAVEGUIDE_WIDTH = 0.5
TAPER_WIDTH = 3
TAPER_LENGTH = 20


# The function below is an reference. You need to provide an EBEAM_TECH
# or replace the layer in the call to layout_waveguide_from_points
def layout_ebeam_waveguide_from_points(
    cell, points_list, radius=None, width=None, taper_width=None, taper_length=None
):
    """ Takes a list of points and lays out a rounded waveguide with optional tapers
    """

    TECHNOLOGY = EBeam
    if radius is None:
        radius = WAVEGUIDE_RADIUS
    if width is None:
        width = WAVEGUIDE_WIDTH
    if taper_width is None:
        taper_width = TAPER_WIDTH
    if taper_length is None:
        taper_length = TAPER_LENGTH

    from zeropdk.layout.waveguide_rounding import layout_waveguide_from_points

    layout_waveguide_from_points(
        cell,
        TECHNOLOGY.layers["Si"],
        points_list,
        width,
        radius,
        taper_width,
        taper_length,
    )

    return cell


def draw_ports(cell, ports):
    """ Draws ports in the Pin Recognition layer (SiEPIC)
    """

    if isinstance(ports, abc.Mapping):  # dictionary
        for port in ports.values():
            port.draw(cell, EBeam.layers["PinRec"])
    elif isinstance(ports, abc.Sequence):  # list
        for port in ports:
            port.draw(cell, EBeam.layers["PinRec"])
    else:
        raise RuntimeError("Give a list or dict of Ports")
