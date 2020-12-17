import pya  # this is the main module of klayout

# system imports
from predefined import MZI_Broadband_DC


def main():
    layout = pya.Layout()
    dbu = layout.dbu = 0.001
    TOP = layout.create_cell("TOP")

    origin = pya.DPoint(0, 0)
    ex = pya.DVector(1, 0)
    ey = pya.DVector(0, 1)

    MZI_Broadband_DC("MZI1x2").place_cell(TOP, origin)

    # from ebeam_pdk import layout_ebeam_waveguide_from_points

    # points_list = [origin, origin + 20 * ex, origin + 20 * ex + 40 * ey]
    # layout_ebeam_waveguide_from_points(TOP, points_list)

    print("Wrote to example_mask.gds")
    layout.write("example_mask.gds")


if __name__ == "__main__":
    main()
