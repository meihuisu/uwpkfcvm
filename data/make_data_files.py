#!/usr/bin/env python

##
#  Builds the data files in the expected format from 
#
# from >>Latitude Longitude  X??    Y??    depth(KM)  Vp/Vs (km/s)
# 
#        35.13359 -120.57782 -64.0  -66.0  0.0        2.60
#
# depth is in increment of 2000m,
#
# Z is depth below sea level. The model is uniformly gridded in X (4 km),
# Y (6 km), and Z (2 km)

import getopt
import sys
import subprocess
import struct
import array
import ssl
import certifi

if sys.version_info.major >= (3) :
  from urllib.request import urlopen
else:
  from urllib2 import urlopen

## at UWPKFCVM/park.xy.latlon (vp) and park.xy.latlon.S (vs)

model = "UWPKFCVM"

dimension_x = 0
dimension_y = 0 
dimension_z = 0

def usage():
    print("\n./make_data_files.py\n\n")
    sys.exit(0)

def download_urlfile(url,fname):
  print("\ndata file:",url,"\n")
  try:
    response = urlopen(url)
    CHUNK = 16 * 1024
    with open(fname, 'wb') as f:
      while True:
        chunk = response.read(CHUNK)
        if not chunk:
          break
        f.write(chunk)
  except:
    e = sys.exc_info()[0]
    print("Exception retrieving and saving model datafiles:",e)
    raise
  return True

def download_urlfile2(url, fname):
    print("\ndata file:", url, "\n")
    try:
        context = ssl.create_default_context(cafile=certifi.where())
        response = urlopen(url, context=context)
        CHUNK = 16 * 1024
        with open(fname, 'wb') as f:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                  break
                f.write(chunk)
    except Exception as e:
        print("Exception retrieving and saving model datafiles:", e)
        raise
    return True


def main():

    # Set our variable defaults.
    path = ""
    mdir = ""

    try:
        fp = open('./config','r')
    except:
        print("ERROR: failed to open config file")
        sys.exit(1)

    ## look for model_data_path and other varaibles
    lines = fp.readlines()
    for line in lines :
        if line[0] == '#' :
          continue
        parts = line.split('=')
        if len(parts) < 2 :
          continue;
        variable=parts[0].strip()
        val=parts[1].strip()

        if (variable == 'model_data_path') :
            path = val + '/' + model
            continue
        if (variable == 'model_dir') :
            mdir = "./"+val
            continue
        if (variable == 'nx') :
            dimension_x = int(val)
            continue
        if (variable == 'ny') :
            dimension_y = int(val)
            continue
        if (variable == 'nz') :
            dimension_z = int(val)
            continue
        continue
    if path == "" :
        print("ERROR: failed to find variables from config file")
        sys.exit(1)

    fp.close()

#    print("\nDownloading model file\n")
#
#    fname="./"+"park.xy.latlon"
#    url = path + "/" + fname
#    download_urlfile2(url,fname)
#    fname="./"+"park.xy.latlon.S"
#    url = path + "/" + fname
#    download_urlfile2(url,fname)

    subprocess.check_call(["mkdir", "-p", mdir])

    # Now we need to go through the data files and put them in the correct
    # format. More specifically, we need a vp.dat

    fvp = open("./park.xy.latlon");
    f_vp = open("./uwpkfcvm/vp.dat", "wb")

    vp_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))

    print ("dimension is", (dimension_x * dimension_y * dimension_z))

    vp_nan_cnt = 0
    vp_total_cnt =0;

    x_pos=0;
    y_pos=0;
    z_pos=0;

    for line in fvp:
        arr = line.split()

        vp = -1.0
        skip_lat = float(arr[0])
        skip_lon = float(arr[1])
        skip_x = float(arr[2])
        skip_y = float(arr[3])
        in_z = float(arr[4])
        tmp = arr[5]

        if( tmp != "NaN" ) :
           vp = float(tmp)
           vp = vp * 1000.0;
        else:
           vp_nan_cnt = vp_nan_cnt + 1

        vp_total_cnt = vp_total_cnt + 1

        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vp_arr[loc] = vp

        x_pos = x_pos + 1
        if(x_pos == dimension_x) :
          x_pos = 0;
          y_pos = y_pos+1
          if(y_pos == dimension_y) :
            y_pos=0;
            z_pos = z_pos+1
            if(z_pos == dimension_z) :
              print ("All DONE")

    vp_arr.tofile(f_vp)

    fvp.close()
    f_vp.close()
    print("Done! with NaN(", vp_nan_cnt, ") total(", vp_total_cnt,")")


    # Now the vs.dat
    fvs = open("./park.xy.latlon.S");
    f_vs = open("./uwpkfcvm/vs.dat", "wb")

    vs_arr = array.array('f', (-1.0,) * (dimension_x * dimension_y * dimension_z))

    print ("dimension is", (dimension_x * dimension_y * dimension_z))

    vs_nan_cnt = 0
    vs_total_cnt =0;

    x_pos=0;
    y_pos=0;
    z_pos=0;

    for line in fvs:
        arr = line.split()

        vs = -1.0
        skip_lat = float(arr[0])
        skip_lon = float(arr[1])
        skip_x = float(arr[2])
        skip_y = float(arr[3])
        in_z = float(arr[4])
        tmp = arr[5]

        if( tmp != "NaN" ) :
           vs = float(tmp)
           vs = vs * 1000.0;
        else:
           vs_nan_cnt = vs_nan_cnt + 1

        vs_total_cnt = vs_total_cnt + 1

        loc =z_pos * (dimension_y * dimension_x) + (y_pos * dimension_x) + x_pos
        vs_arr[loc] = vs

        x_pos = x_pos + 1
        if(x_pos == dimension_x) :
          x_pos = 0;
          y_pos = y_pos+1
          if(y_pos == dimension_y) :
            y_pos=0;
            z_pos = z_pos+1
            if(z_pos == dimension_z) :
              print ("All DONE")

    vs_arr.tofile(f_vs)

    fvs.close()
    f_vs.close()
    print("Done! with NaN(", vs_nan_cnt, ") total(", vs_total_cnt,")")


if __name__ == "__main__":
    main()

