import abfile
import numpy as np
import matplotlib as mpl
# mpl.use('TkAgg')



class Agrid:
    """ Class corresponding to TOPAZ grid objects: regional grids an depth grids"""
    def __init__(self,file,region_name,region_list,idm=800,jdm=760,type='grid',lon=None,lat=None,mask=None):
        """
        Constructor for the Agrid class. The grid is constructed using classes the the file abfile.py, these grids are in the Agrid object  but we also add other elements
        :param file: Name of the input file [str]
        :param region_name: Name of the region of interest [str]
        :param region_list: Lists of available regions and their corresponding extent [dict]
        :param idm: idm [int]
        :param jdm: jdm [int]
        :param type: type of grid; 'grid' for regional file and 'depth' for depth file
        :param lon:
        :param lat:
        :param mask:
        """
        self.__type=type
        self.__region_name=region_name
        self.__Region = region_list


        if type=='depth':
            self.grid = abfile.ABFileBathy(file,'r',idm=idm,jdm=jdm)
            self.depth = self.grid.read_field('depth').flatten()

            self.__idm = idm
            self.__jdm = jdm

            self.lon = lon #In the depth file, we have to use lon/lat from the grid file
            self.lat = lat

            if region_name!='Arctic':
                self.apply_regional_mask(mask=mask)


        elif type=='grid':
            self.grid=abfile.ABFileGrid(file,'r')
            self.lon=self.grid.read_field('plon').flatten()
            self.lat = self.grid.read_field('plat').flatten()

            self.scpx = self.grid.read_field('scpx')
            self.scpy = self.grid.read_field('scpy')

            self.__idm = self.grid.idm
            self.__jdm = self.grid.jdm

            self.TOPAZ_idm=self.grid.idm #Store original values of idm and jdm, useful when applying regional mask
            self.TOPAZ_jdm = self.grid.jdm

            if region_name!='Arctic':
                self.mask=self.apply_regional_mask()
            else: self.mask=None


    def apply_regional_mask(self,mask=None):
        #The regional mask function is only used for debugging

        if self.__type=='grid':
            min_lon,min_lat,max_lon,max_lat = self.__Region[self.__region_name]

            LON,LAT=self.lon.reshape((self.jdm,self.idm)),self.lat.reshape((self.jdm,self.idm))


            mask=np.where((LON<max_lon) & (LON>min_lon) & (LAT<max_lat) & (LAT>min_lat))
            i_min,i_max=np.min(mask[0]),np.max(mask[0])
            j_min, j_max = np.min(mask[1]), np.max(mask[1])


            new_jdm,new_idm=i_max-i_min,j_max-j_min

            self.lon = self.lon.reshape((self.jdm, self.idm))[i_min:i_max, j_min:j_max].flatten()
            self.lat = self.lat.reshape((self.jdm, self.idm))[i_min:i_max, j_min:j_max].flatten()

            self.jdm = new_jdm
            self.idm = new_idm

            # mask = np.where((self.lon < max_lon) & (self.lon > min_lon) & (self.lat < max_lat) & (self.lat > min_lat))
            #
            # self.lon = self.lon[mask]
            # self.lat = self.lat[mask]


            return mask

        if self.__type=='depth':
            i_min,i_max=np.min(mask[0]),np.max(mask[0])
            j_min, j_max = np.min(mask[1]), np.max(mask[1])

            new_jdm,new_idm=i_max-i_min,j_max-j_min

            self.depth=self.depth.reshape((self.jdm, self.idm))[i_min:i_max, j_min:j_max].flatten()

            self.jdm= new_jdm
            self.idm = new_idm

            # self.depth = self.depth[mask] #In the depth file we have to use the mask from the grid file, because we can't make a new one




    @property
    def idm(self):
        """Function to get hidden variable self.__idm"""
        return self.__idm

    @property
    def jdm(self):
        return self.__jdm

    @idm.setter
    def idm(self,i):
        """Function to change variable self.__idm"""
        self.__idm=i

    @jdm.setter
    def jdm(self,j):
        self.__jdm=j