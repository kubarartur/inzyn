import sys
from copy import deepcopy

endnames = ["endx", "endy", "endz"]
inames = ["i", "j", "k"]
b_indices = ["ii", "jj", "kk"]

set_endnames = set(endnames)


def write_head(file):
    file.write("""
submodule(hamiltonian) autogen_vdd_boundary
    implicit none
contains

    subroutine vdd_boundary(arr, gr, dens)
        real(wp), intent(out) :: arr(:,:,:)
        type(grid_params), intent(in) :: gr
        real(wp), intent(in) :: dens(:,:,:)

        real(wp) :: dVoneo4pi
        integer :: {n0}
        integer :: {n1}
        integer :: {n2}
        integer :: {i0}, {i1}, {i2}
        integer :: {b0}, {b1}, {b2}

        dVoneo4pi= 0.25_wp/pi*gr%dr(1)*gr%dr(2)*gr%dr(3)
        """.format(n0=endnames[0], n1=endnames[1], n2=endnames[2],
                   i0=inames[0], i1=inames[1], i2=inames[2],
                   b0=b_indices[0], b1=b_indices[1], b2=b_indices[2]))

    # TODO: Find a better name for obnoxious one4opi
    file.write("\n")

    for i, name in enumerate(endnames):
        file.write(
            " "*8 + "{nm} = size(arr, {ip1}) \n".format(nm=name, ip1=i+1))
    file.write("\n")

    file.write("""
        arr(1,:,:) = 0
        arr(endx,:,:) = 0
        arr(:,1,:) = 0
        arr(:,endy,:) = 0
        arr(:,:,1) = 0
        arr(:,:,endz) = 0
               """)
    file.write("\n")


def write_tail(file):
    file.write("""    end subroutine vdd_boundary
end submodule""")


def write_loop(file, boundary_indices, indent=0):

    closure=""
    not_looped=set()
    for i in range(3):
        if i not in boundary_indices:
            not_looped.add(i)

    # what to write to
    inpts=["arr("]
    arr_inputs=[[]]

    for i in range(3):
        if i in not_looped:
            ini_len=len(inpts)
            inpts.extend(deepcopy(inpts))
            arr_inputs.extend(deepcopy(arr_inputs))
            for j, inp in enumerate(inpts[:ini_len]):
                inpts[j]=inp + "1,"
                arr_inputs[j].append("1")
            for j, inp in enumerate(inpts[ini_len:]):
                inpts[ini_len + j]=inp + "{enm},".format(enm=endnames[i])
                arr_inputs[ini_len + j].append("{enm}".format(enm=endnames[i]))
        else:
            for j, inp in enumerate(inpts):
                inpts[j]=inp + "{bi},".format(bi=b_indices[i])
                arr_inputs[j].append("{bi}".format(bi=b_indices[i]))

    for j, _ in enumerate(inpts):
        inpts[j]=inpts[j][:-1] + ")"

    # loop over the points on boundary
    for index in boundary_indices:
        file.write(
            " "*indent + "do {v1} = 2,{v2} - 1 \n".format(v1=b_indices[index],
                                                          v2=endnames[index])
        )
        closure=closure + " "*indent + "end do\n"
    # loop over the points on the inside
    for i in range(2, -1, -1):
        file.write(
            " "*indent + "do {v1} = 2,{v2} - 1\n".format(v1=inames[i],
                                                         v2=endnames[i]))
        closure=closure + " "*indent + "end do\n"

    # for each grid boundary point in this loops body
    for j, inp in enumerate(inpts):
        # point = point + 1/4/pi/dist(
        file.write(" "*(indent + 4) + inp + " = " +
                   inp + " + dens(")
        for k in range(3):
            if k > 0:
                file.write(" ")

            file.write("{nm}".format(nm=inames[k]))
            if k < 2:
                file.write(",")

        file.write(")*dVoneo4pi*&\n")
        file.write(" "*indent + "gr%inv_dist(")

        # dist indices
        for k in range(3):
            if k > 0:
                file.write(" ")
            if arr_inputs[j][k] == "1":
                file.write("{v1} - {v2}".format(v1=inames[k],
                                                v2=arr_inputs[j][k]))
            elif arr_inputs[j][k] in set_endnames:
                file.write("{v1} - {v2}".format(v1=arr_inputs[j][k],
                                                v2=inames[k]))
            else:
                file.write("abs({v1} - {v2})".format(v1=arr_inputs[j][k],
                                                     v2=inames[k]))
            if k < 2:
                file.write(",")

        file.write(") \n")

    return closure


def main(fname):
    with open(fname, "w") as file:
        write_head(file)

        # all corners
        closure=write_loop(file, [], indent=8)
        file.write(closure)
        file.write("\n")

        # all edges
        for i in range(3):
            # edge parallel to i axis
            closure=write_loop(file, [i], indent=8)
            file.write(closure)
            file.write("\n")

        # all faces
        for i in range(3):
            for j in range(i+1, 3):
                # Faces parallel to X_i X_j plane
                closure=write_loop(file, [j, i], indent=8)
                file.write(closure)
                file.write("\n")

        write_tail(file)


if __name__ == "__main__":
    main(sys.argv[1])
