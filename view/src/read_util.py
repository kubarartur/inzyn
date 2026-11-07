import numpy as np
import f90nml
import src.unit_conv as conv


class Spatial_quantity3D:

    def __init__(self, fname):
        print(f"will read {str(fname) + '.npy'}")
        self.values = np.load(str(fname) + '.npy')
        # verify if the quantity is indeed 3D
        assert (len(self.values.shape) == 3)

        nml = f90nml.read(str(fname) + ".nml")["spacial_params"]
        self.blc = [x/conv.mum_to_atm for x in nml['bottom_left_corner']]
        self.dr = [x/conv.mum_to_atm for x in nml['cell_size']]

    def __str__(self):
        ret = "------------------------------------- \n"
        ret += "spatial 3D quantity with: \n"
        ret += f"shape of values: {self.values.shape} \n"
        ret += f"blc is {self.blc} \n"
        ret += f"dr is {self.dr} \n"

        return ret


class Spatial_quantity2D:

    def __init__(self,
                 quant3D: Spatial_quantity3D,
                 dims, cross=None):

        dims = set(dims)

        slic = []
        self.blc = []
        self.dr = []
        for dim in range(3):
            if dim in dims:
                slic.append(slice(None))
                self.blc.append(quant3D.blc[dim])
                self.dr.append(quant3D.dr[dim])
            else:
                if cross is None:
                    slic.append(quant3D.values.shape[dim]//2)
                    # TODO: mean is the shape is odd?
                else:
                    slic.append(cross)

        slic = tuple(slic)

        self.values = quant3D.values[slic]

    def __str__(self):
        ret = "------------------------------------- \n"
        ret += "spatial 2D quantity with: \n"
        ret += f"shape of values: {self.values.shape} \n"
        ret += f"shape of values: {self.values.shape} \n"
        ret += f"blc is {self.blc} \n"
        ret += f"dr is {self.dr} \n"

        return ret


class Spatial_quantity1D:

    def __init__(self,
                 quant3D: Spatial_quantity3D,
                 dim, crosses=(None, None)):

        slic = []
        self.blc = []
        self.dr = []
        i = 0
        for dim_ in range(3):
            if dim_ == dim:
                slic.append(slice(None))
                self.blc.append(quant3D.blc[dim])
                self.dr.append(quant3D.dr[dim])
            else:
                if crosses[i] is None:
                    slic.append(quant3D.values.shape[dim_]//2)
                    # TODO: mean is the shape is odd?
                else:
                    slic.append(crosses[i])
                i += 1

        slic = tuple(slic)

        self.values = quant3D.values[slic]

    def __str__(self):
        ret = "------------------------------------- \n"
        ret += "spatial 1D quantity with: \n"
        ret += f"shape of values: {self.values.shape} \n"
        ret += f"shape of values: {self.values.shape} \n"
        ret += f"blc is {self.blc} \n"
        ret += f"dr is {self.dr} \n"

        return ret

