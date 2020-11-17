import numpy as np
import chaoticfields as chaos
import matplotlib.pyplot as plt
from shutil import copyfile

def copy_pyd():
    copy_to = 'C:/Users/Max/Documents/Uni/MPhys/chaotic magnetic field/MPhys-Chaotic-fields/chaoticfields.pyd'
    copy_from = 'C:/Users/Max/Documents/Uni/MPhys/chaotic magnetic field/MPhys-Chaotic-Fields-Extension/x64/Release/chaoticfields.pyd'
    copyfile(copy_from, copy_to)
    print('copied')


class Testing:
    def __init__(self):
        self.test_abc()
        self.test_double()
        self.test_wire()

    def test_abc(self):
        l = chaos.abc_field(0, 1, 0.1, [1, 2, 1], [1, 2, 2, 1])
        print(l.s)

    def test_double(self):
        l = chaos.double_abc_field(0, 1, 0.1, [1, 1, 1], [1, 1, 1, 1, 1, 1, 1])

    def test_wire(self):
        l = chaos.wire_field(0, 1, 0.1, [1, 0, 0])

    def test_uniform(self):
        l = chaos.uniform_field(0, 1, 0.1, [1, 0, 0], [1, 0, 0])


if __name__ == '__main__':

    print('nice')
    Testing()
    copy_pyd()

