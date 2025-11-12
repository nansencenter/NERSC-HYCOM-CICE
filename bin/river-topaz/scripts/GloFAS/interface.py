import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets


import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree


class Ui:
    def __init__(self):
        self.app = pg.mkQApp("Parameters")

        params = [
            {'name': 'Input files', 'type': 'group', 'children': [
                {'name': 'TOPAZ grid A file', 'type': 'file', 'value' :'../../data/TOPAZ2/regional.grid.a'},
                {'name': 'TOPAZ depth grid A file', 'type': 'file', 'value' :'../../data/TOPAZ2/regional.depth.a'},
                {'name': 'GloFAS daily disharge file', 'type': 'file', 'value' :'../../data/river/daily_glofas_discharge_files/river_global_06_10.nc'},
                {'name': 'GloFAS LDD file', 'type': 'file', 'value' :'../../data/river/auxiliary_files__glofas_v4_0/ldd_glofas_v4_0.nc'},
                {'name': 'GloFAS uparea file', 'type': 'file', 'value' :'../../data/river/auxiliary_files__glofas_v4_0/uparea_glofas_v4_0.nc'},
                {'name': 'GloFAS elevation file', 'type': 'file', 'value' :'../../data/river/auxiliary_files__glofas_v4_0/elevation_glofas_v4_0.nc'},
                {'name': 'Estuaries edit file', 'type': 'file', 'value': '../../data/river/GloFAS_edit_files/estuaries_edit.csv'},
                {'name': 'River correction file', 'type': 'file', 'value': '../../data/river/correction_files/largest_river_correction_2001_2019.nc'}
            ]},

            {'name': 'Output files', 'type': 'group', 'children': [
                {'name': 'Estuaries file', 'type': 'file', 'value': '../../outputs/GloFAS_estuaries.csv'},
                {'name': 'Modification file', 'type': 'file', 'value': '../../outputs/estuaries_edit.csv'},

            ]},

            {'name': 'Parameters', 'type': 'group', 'children': [
                {'name': 'Region', 'type': 'str', 'value':'Arctic'},
                {'name': 'field name', 'type': 'str', 'value': 'dis24'},
                {'name': 'uparea threshold [m2]', 'type': 'float', 'value': 1e9},
                {'name': 'elevation threshold [m]', 'type': 'float', 'value': 10},
                {'name': 'Month', 'type': 'float', 'value': 1},
                {'name': 'Correct GloFAS data', 'type': 'bool', 'value': True},
                {'name': 'Climatology', 'type': 'bool', 'value': False},
                {'name': 'Edit estuaries', 'type': 'bool', 'value': False}

            ]},
            {'name': 'Run interpolation', 'type': 'group', 'children': [
                {'name': 'Run', 'type': 'action'}]}
        ]

        ## Create tree of Parameter objects
        self.p = Parameter.create(name='params', type='group', children=params)

        self.p.sigTreeStateChanged.connect(self.change)

        self.p.param('Run interpolation', 'Run').sigActivated.connect(self.run)

        ## Create two ParameterTree widgets, both accessing the same data
        t = ParameterTree()
        t.setParameters(self.p, showTop=False)
        t.setWindowTitle('pyqtgraph example: Parameter Tree')

        win = QtWidgets.QWidget()
        layout = QtWidgets.QGridLayout()
        win.setLayout(layout)
        layout.addWidget(QtWidgets.QLabel("Parameters of the model"),
                         0, 0, 1, 2)
        layout.addWidget(t, 1, 0, 1, 1)
        win.show()

        pg.exec()


    ## If anything changes in the tree, print a message
    def change(self,param, changes):
        print("tree changes:")
        for param, change, data in changes:
            path = self.p.childPath(param)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = param.name()
            print('  parameter: %s' % childName)
            print('  change:    %s' % change)
            print('  data:      %s' % str(data))
            print('  ----------')




    # def valueChanging(param, value):
    #     print("Value changing (not finalized): %s %s" % (param, value))

    def run(self):
        self.app.quit()








