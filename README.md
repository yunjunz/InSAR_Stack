# InSAR_Stack

InSAR_Stack is an under-developing set of scripts to prepare for moving PySAR to python3 and with more advanced features. The main new feature so far is the insarStack which is designed to handle networks of InSAR data of different types from different platforms and different orbits, and write the entire dataset to a singe HDF5 file with following structure:

/platform1-track1/ ->	observations (3D) :  unwrapped interferograms, RangeOffset, AzimuthOffset, ...
        		quality  (3D)    :  coherence, uncertainty, ...
        		geometry (2D or 3D)    :  incidence, heading angle, latitute, longitude, ..	

/platform1-track2/ ->   observations
			quality
			geometry

/platform1-track3/ ->   ...

/platform2-track1/ ->   ...

##############################

In a higher level user needs to write scripts to use insarStack.py for writing HDF5 files. In case of one plateform and one single track with multiple datasets for observations, quality and geometry processed with ISCE software, I have included a script "loadData.py". 

%%%%%%%%%%%%%%%%%
NOTE: I have included a simple reader which uses gdal for reading the InSAR data. Therefore we assume that the data are readable by gdal.
%%%%%%%%%%%%%%%%%

Example 1:

loading a stack of unwrapped interferograms to an HDF5 file:

	loadData.py -i /home/hfattahi/Data/T103/interferograms -p filt*unw -d unwrapped-phase -o pairs.h5

loadData.py expects directories of pairs with YYYYMMDD_YYYYMMDD formats inside "/home/hfattahi/Data/T103/interferograms" and looks for name patterns of filt*unw and writes a 3D dataset called "unwrapped-phase" to observations group in pairs.h5 output file.

Note: A 2D dataset (called pairs_idx) which contains index_number of slave and master dates of each pair exists for each platform-track in the output HDF5 file. This way we can always know the dates of each page of the 3D qube dataset.

Example 2:

loading a stack of unwrapped interferograms and estimated ionospheric phase to an HDF5 file:
   	
	loadData.py -i /home/hfattahi/Data/T103/interferograms /home/hfattahi/Data/T103/ionosphere  -p filt*unw iono.bil -d unwrapped-phase iono -o pairs.h5


Example 3:

loading a stack of 
		a) unwrapped interferograms 
		b) estimated ionospheric phase 
		c) coherence of interferograms 
		d) uncertainty of estimated ionospheric phase 
	to an HDF5 file:
        
        loadData.py -i /home/hfattahi/Data/T103/interferograms /home/hfattahi/Data/T103/ionosphere  -p filt*unw iono.bil -d unwrapped-phase iono -I /home/hfattahi/Data/T103/interferograms /home/hfattahi/Data/T103/ionosphere  -P filt*cor Sigma_iono.bil -q  coherence sig_iono -o pairs.h5

Note that the coherence and sig_iono files were considered as quality datasets.


Example 4:

loading a stack of
	observations:
                a) unwrapped interferograms
                b) estimated ionospheric phase
	quality:
                c) coherence of interferograms
                d) uncertainty of estimated ionospheric phase
	geometry:
		e) latitude.rdr
		f) longitude.rdr
		g) height.rdr
		h) waterMask.rdr
		i) los.rdr


	loadData.py -i /home/hfattahi/Data/T103/interferograms /home/hfattahi/Data/T103/ionosphere  -p filt*unw iono.bil -d unwrapped-phase iono -I /home/hfattahi/Data/T103/interferograms /home/hfattahi/Data/T103/ionosphere  -P filt*cor Sigma_iono.bil -q  coherence sig_iono --file_list /home/hfattahi/Data/T103/geometry/lat.rdr /home/hfattahi/Data/T103/geometry/lon.rdr  /home/hfattahi/Data/T103/geometry/z.rdr /home/hfattahi/Data/T103/geometry/waterMask.rdr /home/hfattahi/Data/T103/geometry/los.rdr /home/hfattahi/Data/T103/geometry/los.rdr  --name_list lat lon water-mask height incidence-angle heading-angle -b 1 1 1 1 1 2  -o pairs.h5

Note: with --file_list any 2D dataset can be stored in the HDF5 file. The full path to the exact file should be given. Since the data can be a band in a multi-band data file, then the band number should be specified with option -b




