# =============================================================
# this piece of code was not particularily loved
# =============================================================
import sys


def write_head(file, indices):
    file.write("""
submodule(hamiltonian) autogen_lagrangeian
    implicit none
contains

    subroutine laplaceian(gr, psi_in, psi_out)
        type(grid_params), intent(in) :: gr
        complex(wp), intent(in) :: psi_in(:,:,:)
        complex(wp), intent(out) :: psi_out(:,:,:)

        integer :: endx
        integer :: endy
        integer :: endz
        integer :: """)

    wlist = []
    for i in range(3):
        wlist.append(indices[i])
        wlist.append(", ")

    wlist.pop()
    wlist.append("""
        endx = size(psi_in, 1)
        endy = size(psi_in, 2)
        endz = size(psi_in, 3)
""")
    wlist.append("\n")
    wlist.append("\n")

    for ws in wlist:
        file.write(ws)


def write_tail(file):
    file.write("""
    end subroutine laplaceian
end submodule""")


def write_psi_input(file, indices, indent=0):
    psilist = [" "*indent + "psi_out("]

    for indx in indices:
        psilist.extend(indx)
        psilist.append(",")

    # pop the trailing comma
    psilist.pop()

    psilist.append(")")
    psilist.append(" = &\n")

    for stri in psilist:
        file.write(stri)


def write_loop_head(file, i, iname, indent=0):
    wstr = " "*indent + "do " + iname + \
        " = 2,size(psi_in, " + str(i + 1) + ") - 1\n"
    file.write(wstr)


def write_central_diff(file, inames, der_dir, indent=0):
    wstr = [" "*indent + "("]

    # first term
    wstr.append("psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(name)
            wstr.append("+1")
            wstr.append(",")
        else:
            wstr.append(name)
            wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # second term
    wstr.append(" - ")
    wstr.append("2._wp*psi_in(")
    for i, name in enumerate(inames):
        wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # third term
    wstr.append(" + ")
    wstr.append("psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(name)
            wstr.append("-1")
            wstr.append(",")
        else:
            wstr.append(name)
            wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # closing division
    wstr.append(")/gr%dr(")

    # due to fortran indexing issues
    wstr.append(str(der_dir + 1))
    wstr.append(")**2")

    for ws in wstr:
        file.write(ws)


def write_forward_diff(file, inames, der_dir, indent=0):
    wstr = [" "*indent + "("]

    i_start = int(inames[der_dir])

    # 1/4 term
    wstr.append("2._wp*psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(str(name))
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 2/4 term
    wstr.append(" - ")
    wstr.append("5._wp*psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(str(i_start + 1))
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 3/4 term
    wstr.append(" + &\n")
    wstr.append(" "*indent)
    wstr.append("4._wp*psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(str(i_start + 2))
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 4/4 term
    wstr.append(" - ")
    wstr.append("psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(str(i_start + 3))
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # closing division
    wstr.append(")/gr%dr(")

    # due to fortan indexing issues
    wstr.append(str(der_dir + 1))
    wstr.append(")**2")

    for ws in wstr:
        file.write(ws)


def write_backward_diff(file, inames, der_dir, indent=0):
    wstr = [" "*indent + "("]

    i_start = inames[der_dir]

    # 1/4 term
    wstr.append("2._wp*psi_in(")
    for i, name in enumerate(inames):
        wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 2/4 term
    wstr.append(" - ")
    wstr.append("5._wp*psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(i_start + "-1")
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 3/4 term
    wstr.append(" + &\n")
    wstr.append(" "*indent)
    wstr.append("4._wp*psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(i_start + "-2")
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # 4/4 term
    wstr.append(" - ")
    wstr.append("psi_in(")
    for i, name in enumerate(inames):
        if i == der_dir:
            wstr.append(i_start + "-3")
        else:
            wstr.append(name)
        wstr.append(",")
    wstr.pop()
    wstr.append(")")

    # closing division
    wstr.append(")/gr%dr(")

    # due to fortan indexing issues
    wstr.append(str(der_dir + 1))
    wstr.append(")**2")

    for ws in wstr:
        file.write(ws)


def write_inner_outer(file, indices, signs):
    indices_inner = ["nan", "nan", "nan"]
    end_names = ["endx", "endy", "endz"]

    for i, name in enumerate(indices):
        if (signs[i] == "i"):
            indices_inner[i] = indices[i]
        elif (signs[i] == "-"):
            indices_inner[i] = end_names[i]
        elif (signs[i] == "+"):
            # TODO: this should not break when chaning from string to integer
            indices_inner[i] = "1"

    for i, name in enumerate(indices):
        if (signs[i] == "i"):
            write_loop_head(file, i, indices[i], indent=8)

    write_psi_input(file, indices_inner, indent=12)
    for i in range(3):
        if signs[i] == "i":
            write_central_diff(file, indices_inner, i, indent=16)
        elif signs[i] == "+":
            write_forward_diff(file, indices_inner, i, indent=16)
        elif signs[i] == "-":
            write_backward_diff(file, indices_inner, i, indent=16)
        else:
            raise ValueError(f"Expected +,- or i as sign got : {signs[i]} ")

        if i < 2:
            file.write("+ &\n")

    file.write("\n")
    for i in range(3):
        if signs[i] == "i":
            file.write(" "*8 + "end do\n")


def main(fout):
    with open(fout, "w") as file:
        indices = ["ii", "jj", "kk"]
        write_head(file, indices)
        # actuall logic here
        signs_arr = [
            ["i", "i", "i"],

            ["+", "i", "i"],
            ["-", "i", "i"],

            ["i", "+", "i"],
            ["i", "-", "i"],

            ["i", "i", "+"],
            ["i", "i", "-"],

            ["+", "+", "i"],
            ["+", "-", "i"],

            ["-", "+", "i"],
            ["-", "-", "i"],

            ["+", "+", "+"],
            ["+", "+", "-"],
            ["+", "-", "+"],
            ["+", "-", "-"],

            ["-", "+", "+"],
            ["-", "+", "-"],
            ["-", "-", "+"],
            ["-", "-", "-"],
        ]
        for signs in signs_arr:
            write_inner_outer(file, indices, signs)
            file.write("\n")

        write_tail(file)


if __name__ == "__main__":
    main(sys.argv[1])
