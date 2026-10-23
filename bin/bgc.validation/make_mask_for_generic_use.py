import numpy as np

def get_mask(x,y):
  # if your lat and lon are 1D, do the following before calling get_mask
  #x,y = np.meshgrid(lon,lat)

  maskall = np.zeros((y.shape))

  # barents (BRTS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>16.), x<=58.)] = 1
  brts=masks.copy()
  maskall = maskall + masks*10.0

  # lofoten (LFTB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>7.), x<=16.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<73.), x>-10.), x<=7.)] = 1
  lftb=masks.copy()
  maskall = maskall + masks*1.0

  # norwegian (NRWB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=60 , y<68.), x>-10.), x<=16.)] = 1
  nrwb=masks.copy()
  maskall = maskall + masks*6.0

  # irminger (IBIB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=55 , y<66.), x>-46.), x<=-10.)] = 1
  ibib=masks.copy()
  maskall = maskall + masks*3.0

  # greenland (GSIS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=73 , y<81.), x>-25.), x<=7.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=66 , y<73.), x>-37.), x<=-10.)] = 1
  gsis=masks.copy()
  maskall = maskall + masks*13.0

  # kara (KARA)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>58.), x<=100.)] = 1
  kara=masks.copy()
  maskall = maskall + masks*12.0

  # eurasia (ERSB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=78 , y<90.), x>100.), x<=140.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=81 , y<90.), x>-25.), x<=100.)] = 1
  ersb=masks.copy()
  maskall = maskall + masks*2.0

  # laptev (LPES)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<78.), x>100.), x<=180.)] = 1
  lpes=masks.copy()
  maskall = maskall + masks*8.0

  # amerasia (ARSB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=81 , y<90.), x>-123.), x<=-25.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=78 , y<90.), x>140.), x<=180.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=76 , y<90.), x>-181.), x<=-123.)] = 1

  arsb=masks.copy()
  maskall = maskall + masks*14.0

  # chukchi (CHBF)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=66 , y<76.), x>-180.), x<=-158.)] = 1
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<76.), x>-158.), x<=-123.)] = 1

  chbf=masks.copy()
  maskall = maskall + masks*4.0

  # bering (BRGS)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=56 , y<66.), x>=-180.), x<=-158.)] = 1

  brgs=masks.copy()
  maskall = maskall + masks*9.0

  # canadian (CACH)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>-123.), x<=-78.)] = 1

  cach=masks.copy()
  maskall = maskall + masks*5.0

  # baffin (BFFB)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=68 , y<81.), x>-78.), x<=-51.)] = 1

  bffb=masks.copy()
  maskall = maskall + masks*11.0

  # labrador (LBRD)
  masks = np.zeros((y.shape))
  masks[np.logical_and(np.logical_and(np.logical_and(y>=55 , y<68.), x>-64.), x<=-46.)] = 1

  lbrd=masks.copy()
  maskall = maskall + masks*7.0


  masks = {'LBRD':lbrd, 'BFFB':bffb, 'CACH':cach,\
         'BRGS':brgs, 'CHBF':chbf, 'ARSB':arsb,\
          'LPES':lpes, 'ERSB':ersb, 'KARA':kara,\
          'GSIS':gsis, 'IBIB':ibib, 'NRWB':nrwb,\
          'LFTB':lftb, 'BRTS':brts} 

  return masks
