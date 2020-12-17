from ebeam_pdk import PDKCell, GDSCell
from zeropdk.layout.waveguides import layout_waveguide
from zeropdk.layout import insert_shape
from zeropdk.pcell import Port
from zeropdk.pcell import place_cell, ParamContainer
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


def define_param(name, type, description, default=None, **kwargs):
    from zeropdk.pcell import PCellParameter

    return PCellParameter(
        name=name, type=type, description=description, default=default, **kwargs
    )


class Broadband_DC_te1550(GDSCell("ebeam_bdc_te1550", filename="ebeam_bdc_te1550.gds")):
    """
    The PCell version of directional coupler
    """

    def draw(self, cell):
        cell = self.draw_gds_cell(cell)
        _, ex, ey = self.origin_ex_ey()

        opt1_position = -35.45 * ex + 2.35 * ey
        opt2_position = -35.45 * ex - 2.35 * ey
        opt3_position = 35.3 * ex + 2.35 * ey
        opt4_position = 35.3 * ex - 2.35 * ey

        ports = [
            Port("opt1", opt1_position, -ex, 0.5),
            Port("opt2", opt2_position, -ex, 0.5),
            Port("opt3", opt3_position, ex, 0.5),
            Port("opt4", opt4_position, ex, 0.5),
        ]

        return cell, {port.name: port for port in ports}


class YBranch_te1550(GDSCell("ebeam_y_1550", filename="ebeam_y_1550.gds")):
    """
    The PCell version of YBranch, mainly for use in scripts.
    """

    def draw(self, cell):
        cell = self.draw_gds_cell(cell)
        _, ex, ey = self.origin_ex_ey()

        opt1_position = -7.4 * ex + 0 * ey
        opt2_position = 7.4 * ex + 2.75 * ey
        opt3_position = 7.4 * ex + -2.75 * ey

        ports = [
            Port("opt1", opt1_position, -ex, 0.5),
            Port("opt2", opt2_position, ex, 0.5),
            Port("opt3", opt3_position, ex, 0.5),
        ]

        return cell, {port.name: port for port in ports}


class Waveguide_heater(PDKCell):
    params = ParamContainer(
        define_param(
            "mh_length", TypeInt, "Heated waveguide length", default=100, unit="um"
        ),
        define_param("mh_width", TypeDouble, "Heater width", default=5, unit="um"),
        define_param("w_contacts", TypeDouble, "Contact width", default=10, unit="um"),
    )

    def draw(self, cell):
        layout = cell.layout()

        cp = self.get_cell_params()
        lay = self.get_layers(layout)
        origin, ex, ey = self.origin_ex_ey()

        w_mh = cp.mh_width
        length = cp.mh_length
        w_contacts = cp.w_contacts

        input_port_position = origin - length / 2 * ex
        output_port_postion = origin + length / 2 * ex
        x_off = length / 2 - w_contacts / 2
        y_off = w_mh / 2 + w_contacts / 2

        input_contact_position = origin - x_off * ex + y_off * ey
        output_contact_position = origin + x_off * ex + y_off * ey

        # draw metal heater
        layout_waveguide(
            cell, lay.M_Heater, [input_port_position, output_port_postion], w_mh
        )

        # draw metal contacts

        # Left contacts
        contact_points = [
            origin - x_off * ex + (y_off - w_contacts / 2) * ey,
            origin - x_off * ex + (y_off + w_contacts / 2) * ey,
        ]

        layout_waveguide(cell, lay.M_Heater, contact_points, w_contacts)
        layout_waveguide(cell, lay.ML, contact_points, w_contacts)

        # Right contacts
        contact_points = [
            origin + x_off * ex + (y_off - w_contacts / 2) * ey,
            origin + x_off * ex + (y_off + w_contacts / 2) * ey,
        ]
        layout_waveguide(cell, lay.M_Heater, contact_points, w_contacts)
        layout_waveguide(cell, lay.ML, contact_points, w_contacts)

        ports = [
            Port("el_contact_1", input_contact_position, ey, w_contacts),
            Port("el_contact_2", output_contact_position, ey, w_contacts),
        ]

        return cell, {port.name: port for port in ports}


class MZI_Broadband_DC(Broadband_DC_te1550, YBranch_te1550, Waveguide_heater):
    params = ParamContainer(
        define_param(
            "MZI_height", TypeDouble, "interferometer_height of MZI", default=20
        ),
        define_param(
            "waveguide_length_MZI", TypeDouble, "waveguide_length_MZI", default=100
        ),
        define_param("layout_ports", TypeBoolean, "Layout Pins?", default=True),
        define_param("wg_width", TypeInt, "Waveguide width ", default=0.5),
    )

    def draw(self, cell):
        layout = cell.layout()

        cp = self.get_cell_params()
        lay = self.get_layers(layout)
        origin, ex, ey = self.origin_ex_ey()

        YBranch_te1550_Cell_0, YBranch_te1550_Ports_0 = YBranch_te1550(
            "Broadband_DC", params={"angle_ex": cp.angle_ex}
        ).new_cell(layout)
        YBranch_te1550_Ports_0 = place_cell(
            cell, YBranch_te1550_Cell_0, YBranch_te1550_Ports_0, origin
        )

        Broadband_DC_Cell, Broadband_DC_Ports = Broadband_DC_te1550(
            "Broadband_DC", cp
        ).new_cell(layout)
        Broadband_DC_Ports = place_cell(
            cell,
            Broadband_DC_Cell,
            Broadband_DC_Ports,
            origin + (cp.waveguide_length_MZI + 100) * ex,
        )

        top_arm_left = (
            YBranch_te1550_Ports_0["opt2"].position + 25 * ex + cp.MZI_height / 2 * ey
        )
        top_arm_right = top_arm_left + cp.waveguide_length_MZI * ex
        layout_waveguide(cell, lay.Si, [top_arm_left, top_arm_right], 0.5)

        Heater_1_MZI_Cell, Heater_1_MZI_Ports = Waveguide_heater(
            "waveguide heater",
            params={"mh_length": cp.waveguide_length_MZI, "angle_ex": cp.angle_ex},
        ).new_cell(layout)
        Heater_1_MZI_Ports = place_cell(
            cell,
            Heater_1_MZI_Cell,
            Heater_1_MZI_Ports,
            YBranch_te1550_Ports_0["opt2"].position
            + 25 * ex
            + cp.MZI_height / 2 * ey
            + cp.waveguide_length_MZI / 2 * ex,
        )

        bottom_arm_left = (
            YBranch_te1550_Ports_0["opt3"].position + 25 * ex - cp.MZI_height / 2 * ey
        )
        bottom_arm_right = bottom_arm_left + cp.waveguide_length_MZI * ex
        layout_waveguide(cell, lay.Si, [bottom_arm_left, bottom_arm_right], 0.5)

        from zeropdk.layout.geometry import bezier_optimal

        def layout_mzi_curve(P0, P1):
            curve = bezier_optimal(P0, P1, cp.angle_ex, cp.angle_ex)
            return layout_waveguide(cell, lay.Si, curve, cp.wg_width, smooth=True)

        layout_mzi_curve(YBranch_te1550_Ports_0["opt2"].position, top_arm_left)
        layout_mzi_curve(YBranch_te1550_Ports_0["opt3"].position, bottom_arm_left)

        layout_mzi_curve(top_arm_right, Broadband_DC_Ports["opt1"].position)
        layout_mzi_curve(bottom_arm_right, Broadband_DC_Ports["opt2"].position)

        ports = []

        ports.append(YBranch_te1550_Ports_0["opt1"].rename("opt_in"))
        ports.append(Broadband_DC_Ports["opt3"].rename("opt_out_1"))
        ports.append(Broadband_DC_Ports["opt4"].rename("opt_out_2"))
        ports.extend(Heater_1_MZI_Ports.values())

        if cp.layout_ports:
            for port in ports:
                port.draw(cell, lay.PinRec)
            insert_shape(cell, lay.DevRec, cell.bbox())

        return cell, {port.name: port for port in ports}
