# ########################################
# Script to correct continuum ALMA images
# for primary beam. 
##########################################
continuum_dir = 'continuum'
sources = ['16263_2422','16267_2417','AN6','B2-A7','SM1','SM1N']
file_ext = ['_selfcal_2','_taper','_taper','_taper','_selfcal_3','']

for i in range(len(sources)):
    	filename = '{0}/{1}_continuum{2}.image'.format(continuum_dir,sources[i],file_ext[i])
	outfile = '{0}/{1}_continuum{2}.pbcor'.format(continuum_dir,sources[i],file_ext[i])
        pbimage = '{0}/{1}_continuum{2}.flux'.format(continuum_dir,sources[i],file_ext[i])
	impbcor(imagename=filename,
                pbimage=pbimage,
		outfile=outfile)
	exportfits(imagename=outfile,
                   fitsimage='{0}.fits'.format(outfile))



