#!/usr/bin/env python
import argparse
import hashlib
import logging
import tarfile
import shutil
import gzip
import os

# Set up logger
_loglevel=logging.DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(_loglevel)
formatter = logging.Formatter("%(asctime)s - %(name)10s - %(levelname)7s: %(message)s")
ch = logging.StreamHandler()
ch.setLevel(_loglevel)
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate=False

_magic_dict = {
   "gz"  : "\x1f\x8b\x08",
   "bz2" : "\x42\x5a\x68",
   "zip" : "\x50\x4b\x03\x04"
    }
_max_len = max(len(x) for x in _magic_dict)


class TarChunk(object) :

   def __init__(self,prefix,chunk_size,unitfac,unit) :
      self._prefix=prefix
      self._chunk_counter=1
      self._chunk_size=chunk_size
      self._fnametmpl = "%s-%05d.tar"
      self._current_size=0
      self._unitfac=unitfac
      self._unit=unit
      self._tf=None
      self._md5f=None


   def new_filename(self) :
      fname = "%s-%05d.tar"%(self._prefix,self._chunk_counter)
      while os.path.exists(fname) :
         self._chunk_counter=self._=self._chunk_counter+1
         fname = "%s-%05d.tar"%(self._prefix,self._chunk_counter)
      return fname

   def add_file(self,infile,file_size) :
      #logger.info("%s - Size %14.2f%s. Adding %s"%(tf.name,float(current_size)/unitfac,unit,infile))
      logger.info("%s - Size %14.2f%s. Adding %s"%(os.path.basename(self._tf.name),float(self._current_size)/self._unitfac,self._unit,infile))
      self._tf.add(infile,recursive=False)

      ## md5 sum of file
      mymd5sum = hashlib.md5(open(infile, 'rb').read()).hexdigest()
      self._md5f.write("%s  %s\n"%(infile, mymd5sum))
      self._current_size=self._current_size+file_size


   def next_chunk(self) :

      # Close old file
      if self._tf is not None : 
         self._tf.close() 
         self._md5f.close()

      # Next filename
      fname = self.new_filename()
      
      # Open file and set info
      self._tf=tarfile.TarFile(fname,"w")
      self._md5f=open(fname+".md5","w")
      self._current_size=0

   def close() : 
      self._tf.close()
      self._tf=None
      self._md5f=None


   @property
   def chunk_size(self) : return self._chunk_size

   @property
   def current_size(self) : return self._current_size


         


#def new_chunk(chunk_counter,prefix) :
#
#   # Dont overwrite local files
#   fname = "%s-%05d.tar"%(prefix,chunk_counter)
#   while os.path.exists(fname) :
#      fname = "%s-%05d.tar"%(prefix,chunk_counter)
#      chunk_counter=chunk_counter+1
#   
#   # Open file and set info
#   tf=tarfile.TarFile("%s-%05d.tar"%(prefix,chunk_counter),"w")
#   chunk_counter = chunk_counter + 1
#   current_size=0
#   return tf,chunk_counter,current_size
#
#def add_file(tf,infile,current_size,file_size,unitfac,unit) :
#   #logger.info("%s - Size %14.2f%s. Adding %s"%(tf.name,float(current_size)/unitfac,unit,infile))
#   logger.info("%s - Size %14.2f%s. Adding %s"%(os.path.basename(tf.name),float(current_size)/unitfac,unit,infile))
#   tf.add(infile,recursive=False)
#   current_size=current_size+file_size
#   return current_size
#
## Not 100% bomb-proof, but should normally be good enough. The worst that can happen
## is that we avoid compressing a file by mistake.
#def file_is_gzip(fname) :
#   try :
#      tst=open(fname).read(_max_len)
#      if tst.startswith(_magic_dict["gz"]) :
#         return True
#      else :
#         return False
#   except :
#      return False


def get_chunk_size(chunk_size) :
   if chunk_size is None :
      return 1024**3,"G",None

   # Try if this is an integer - interpreted as bytes
   try :
      chunk_size=int(chunk_size)
      fac=1
      unit="B"
   except :
      unit = chunk_size[-1] 
      if unit == "K" :
         fac=1024
      elif unit == "M" :
         fac=1024*1024
      elif unit == "G" :
         fac=1024*1024*1024
      elif unit == "T" :
         fac=1024**4
      else :
         msg="Unable to interpret chunk_size=%s"%chunk_size
         logger.error(msg)
         raise ValueError,msg
      try :
         chunk_size = int(chunk_size[:-1])
         chunk_size=chunk_size * fac
      except :
         msg="Unable to interpret chunk_size=%s"%chunk_size
         logger.error(msg)
         raise ValueError,msg
   return fac,unit,chunk_size



def main(files,prefix="ARCHIVE",compression="",chunk_size=None) :


   unitfac,unit,chunk_size=get_chunk_size(chunk_size)
   logger.info("Chunk_size=%d"%chunk_size)

   tc = TarChunk(prefix,chunk_size,unitfac,unit)
   tc.next_chunk()


   for file in files :

      # Do compression on file-by-file basis. This makes it faster to retrieve/list
      # individual files from a large archive than if the whole tar file is compressed...
      if not compression  :
         infile=file
      elif compression == "gzip" :
         if file_is_gzip(file) :
            infile=file
            logger.info("Not compressing file %s (already gzipped)"%infile)
         else :
            infile="%s.gz"%file
            logger.info("Compressing %s to %s"%(file,infile))
            with open(file, 'rb') as f_in, gzip.open(infile, 'wb') as f_out:
               shutil.copyfileobj(f_in, f_out)
      else :
         raise NotImplementedError,""

      ## md5 sum of file
      #mymd5sum = hashlib.md5(open(infile, 'rb').read()).hexdigest()
      #logger.info("MD5 sum of %s:%s"%(infile,mymd5sum))
      #md5sums.append(mymd5)

      # Get size of next file to be archived
      file_size = os.stat(infile).st_size

      # 1) File fits within tar size. Write as normal
      if tc.chunk_size is None or tc.current_size + file_size < tc.chunk_size :
         current_size = tc.add_file(infile,file_size)
      # 2) File does not fit within tar size
      else :

         # If tarfile has been written to, close existing, open a new one and write to it
         if current_size > 0 :
            tc.next_chunk()
            tc.add_file(infile,file_size)

         # This implies that file_size > chunk size. We must write somewhere, so write to 
         # existing file, then open a new file
         else :
            tc.add_file(infile,file_size)
            tc.next_chunk()


   # Close last file. TODO: Make sure its not empty
   tc.close()



if __name__ == "__main__"  :


   parser = argparse.ArgumentParser(description='')
   parser.add_argument('--chunk_size',        type=str,default=None)
   parser.add_argument('--prefix',      action="store_true", default="ARCHIVE")
   parser.add_argument('--compression', type=str, default="")
   parser.add_argument('filename',   help="",nargs="+")
   args = parser.parse_args()


   main(args.filename,prefix=args.prefix,compression=args.compression,chunk_size=args.chunk_size)
