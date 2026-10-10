import sys
from copy import deepcopy

endnames = ["endx", "endy", "endz"]
inames = ["i", "j", "k"]
b_indices = ["ii", "jj", "kk"]

set_endnames = set(endnames)


def write_head(file):
    file.write(f"""
submodule(hamiltonian) autogen_vdd_boundary
    implicit none
contains

    subroutine vdd_boundary(arr, gr, dens)
        real(wp), intent(out) :: arr(:,:,:)
        type(grid_params), intent(in) :: gr
        real(wp), intent(in) :: dens(:,:,:)

        integer :: {endnames[0]}
        integer :: {endnames[1]}
        integer :: {endnames[2]}
        integer :: {inames[0]}, {inames[1]}, {inames[2]}
        integer :: {b_indices[0]}, {b_indices[1]}, {b_indices[2]}
        """)

    file.write("\n")

    for i, name in enumerate(endnames):
        file.write(" "*8 + f"{name} = size(arr, {i+1}) \n")
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

    closure = ""
    not_looped = set()
    for i in range(3):
        if i not in boundary_indices:
            not_looped.add(i)

    # what to write to
    inpts = ["arr("]
    arr_inputs = [[]]

    for i in range(3):
        if i in not_looped:
            ini_len = len(inpts)
            inpts.extend(deepcopy(inpts))
            arr_inputs.extend(deepcopy(arr_inputs))
            for j, inp in enumerate(inpts[:ini_len]):
                inpts[j] = inp + "1,"
                arr_inputs[j].append("1")
            for j, inp in enumerate(inpts[ini_len:]):
                inpts[ini_len + j] = inp + f"{endnames[i]},"
                arr_inputs[ini_len + j].append(f"{endnames[i]}")
        else:
            for j, inp in enumerate(inpts):
                inpts[j] = inp + f"{b_indices[i]},"
                arr_inputs[j].append(f"{b_indices[i]}")

    for j, _ in enumerate(inpts):
        inpts[j] = inpts[j][:-1] + ")"

    # loop over the points on boundary
    for index in boundary_indices:
        file.write(
            " "*indent + f"do {b_indices[index]} = 2,{endnames[index]} \n")
        closure = closure + " "*indent + "end do\n"

    # loop over the points on the inside
    for i in range(2, -1, -1):
        file.write(" "*indent + f"do {inames[i]} = 2,{endnames[i]} - 2\n")
        closure = closure + " "*indent + "end do\n"

    # for each grid boundary point in this loops body
    for j, inp in enumerate(inpts):
        # point = point + 1/4/pi/dist(
        file.write(" "*(indent + 4) + inp + " = " +
                   inp + " + dens(")
        for k in range(3):
            if k > 0:
                file.write(" ")

            file.write(f"{inames[k]}")
            if k < 2:
                file.write(",")

        file.write(")*0.25_wp/pi/&\n")
        file.write(" "*indent + "gr%dist(")

        # dist indices
        for k in range(3):
            if k > 0:
                file.write(" ")
            if arr_inputs[j][k] == "1":
                file.write(f"{inames[k]} - {arr_inputs[j][k]}")
            elif arr_inputs[i][k] in set_endnames:
                file.write(f"{arr_inputs[j][k]} - {inames[k]}")
            else:
                file.write(f"abs({arr_inputs[j][k]} - {inames[k]})")
            if k < 2:
                file.write(",")

        file.write(") \n")

    return closure


def main(fname):
    with open(fname, "w") as file:
        write_head(file)

        # all corners
        closure = write_loop(file, [], indent=8)
        file.write(closure)
        file.write("\n")

        # all edges
        for i in range(3):
            # edge parallel to i axis
            closure = write_loop(file, [i], indent=8)
            file.write(closure)
            file.write("\n")

        # all faces
        for i in range(3):
            for j in range(i+1, 3):
                # Faces parallel to X_i X_j plane
                closure = write_loop(file, [j, i], indent=8)
                file.write(closure)
                file.write("\n")

        write_tail(file)


if __name__ == "__main__":
    main(sys.argv[1])
