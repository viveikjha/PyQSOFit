'''
This small sample run for the code PyQSOFit basically demonstrates the fitting for a quasar spectra by using QSOFit. Will add a lot of modifications to it in the future depending on our needs!!!!

We can get the fitting parameters, as well as the output physical quantities we obtain from the fit.
'''

import glob, os,sys,timeit,glob
import matplotlib
import numpy as np
from PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings

path= os.getcwd()
path4=path3=path2=path1=path
#This will import  the spectrum from anywhere we want...obviously can load plenty of files (glob)!!!

files=sorted(glob.glob("main_spectra/spec*fits"))
dummy=np.zeros(len(files))

#This will do the fitting one by one.
with open('fitting_results_tmp.txt', 'w') as f:
    for nfiles in range (0,len(files)):
        data=fits.open(files[nfiles])
        lam=10**data[1].data['loglam']        # observed wavelength
        flux=data[1].data['flux']             # observed fluc
        err=1./np.sqrt(data[1].data['ivar'])  # 1 sigma error
        z=data[2].data['z'][0]  # redshift
        ra=data[2].data['PLUG_RA'][0] # RA
        dec=data[2].data['PLUG_DEC'][0] # DEC
        plateid = data[0].header['plateid']   # SDSS plate ID
        mjd = data[0].header['mjd']           # SDSS MJD
        fiberid = data[0].header['fiberid']   # SDSS fiber ID
        q = QSOFit(lam, flux, err, z)
        #print("plate mjd fiber",plateid,mjd,fiberid)
        #print("RA dec",ra,dec)
        start = timeit.default_timer()

    # do the fitting

        q.Fit(name = None,nsmooth = 1, and_or_mask = False, deredden = True, reject_badpix = False, wave_range = None,
        wave_mask =None, decomposition_host =True, Mi = None, npca_gal = 10, npca_qso = 20,
        Fe_uv_op = True, poly = True, BC = False, rej_abs = False, initial_guess = None, MC = False,
        n_trails = 0, linefit = True, tie_lambda = True, tie_width = True, tie_flux_1 = True, tie_flux_2 = True,
        save_result = True, plot_fig = True,save_fig = False, plot_line_name = True, plot_legend = True,
        dustmap_path = path4, save_fig_path =path+'/fit_results/fitted_model_spec-{}-{}-{}'.format(plateid,mjd,fiberid), save_fits_path = path,save_fits_name = 'fit_results/fitting_result_spec-{}-{}-{}'.format(plateid,mjd,fiberid))

        end = timeit.default_timer()
        #print ('Fitting finished in : '+str(np.round(end-start))+'s')




        # To print the details of the fitted parameters.

    # print(q.line_result_name[12:15])
        #print(q.line_result[12:15])
        # Get the emission line parameters (narrow)  { FWHM, area, EW, line centre}
        fwhmo,sigmao,ewo,peako,areao = q.line_prop(4862.68,q.gauss_result[12:15],'narrow')
        fwhmh,sigmah,ewh,peakh,areah = q.line_prop(4862.68,q.gauss_result[0:39],'narrow')
        #print(fwhmo,areao,peako)
        #print('H-beta narrow components',fwhmh,areah,peakh)
        try:
            ratio_n=np.log10(areao/areah) # narrow band flux ratio OIII/H-beta

        except ZeroDivisionError:
            ratio_n=0
            pass



        #Get the emission line parameters (broad) { FWHM, area, EW, line centre}
        t=float(q.line_result[14:15]) # broad fwhm
        s=float(q.line_result[17:18]) # very broad fwhm
        aha=q.line_result[66:67] # area of h-alpha
        fwhm,sigma,ew,peak,area = q.line_prop(q.linelist[6][0],q.gauss_result[0:27],'broad')
        print(q.CalFWHM(8.35,t),q.CalFWHM(8.35,s),fwhm,fwhmh)
        if (area==0):
            area=0.01
        try:
            ratio_b=float(aha)/float(area)

        except TypeError:
            ratio_b=0
            pass


        # Luminosity calculation

        blr_size=10**(1.527+(0.533*np.log10((10**q.conti_result[16])/(10**44))))
         # blr size from the R-L relation of bentz.




    # print('BLR size',blr)
    # print('SMBH mass',mbh)

    #plot the quasar rest frame spectrum after removed the host galaxy component
        #
        #plt.plot(q.wave,q.err,'r')
        #plt.plot(q.wave,q.line_flux,'m',lw=2)
        #plt.plot(q.wave,q.line_fit,'c',lw=2)

    #To plot the whole model, we use Manygauss to reappear the line fitting results saved in gauss_result


        #
        #plt.plot(q.wave,q.Manygauss(np.log(q.wave),q.gauss_result)+q.f_conti_model,'b')
    # plt.plot(q.wave,q.Manygauss(np.log(q.wave),q.gauss_result))
        #print(q.wave,q.Manygauss(np.log(q.wave),q.gauss_result))



        fe_blend=q.f_conti_model-q.PL_poly_BC
        fe_b=np.where((q.wave > 4433) & (q.wave < 4684))
        #plt.plot(q.wave[fe_b],fe_blend[fe_b])
        integrated=np.trapz(fe_blend[fe_b],q.wave[fe_b]) # rfe calculation
        '''
        print('SNR:',q.SN_ratio_conti)
        print('H_beta broad component FWHM:',q.CalFWHM(0,t))
        print('H_beta very broad FWHM:',q.CalFWHM(0,s))
        print('Equivalend width H_beta:',ew)
        print('Area H_b:',area)
        print('OIII/ H_b (narrow component) ratio:',ratio_n)
        print('H_a/ H_b (broad component) ratio:',ratio_b)
        print('FeII to H-beta ratio (R_fe)',(integrated/area))
        print(q.conti_result_name[16],q.conti_result[16])
        print('BLR size (light days):',blr_size)
        '''

        # To calculate the asymmetry in H-beta profile.

        print('no %i done'%nfiles)

        hbeta_lam=np.where((q.wave > 4780) & ( q.wave< 4950))
        ss=q.Manygauss(np.log(q.wave),q.gauss_result)
        hbeta_flux=ss[hbeta_lam]
        hbeta_range=q.wave[hbeta_lam]
        #print(len(hbeta_range),len(hbeta_flux))
        max_flux=np.max(hbeta_flux)
        eig=0.75*max_flux
        twe=0.25*max_flux

        #k=np.where (hbeta_flux ==max_flux)
        eg=np.where((hbeta_flux > 0.9*eig) & (hbeta_flux < 1.1*eig))
        temp2=hbeta_range[eg]
        if(len(temp2) < 2):
                    print("asymmetry can't be calculated at 80%")
                    continue

        for i in range(len(temp2)):

                if(len(temp2)==2):
                    lam_hb,lam_hr=(temp2[i],temp2[i-1])
                else:
                    temp3= temp2[i]-temp2[i-1]
                    if (temp3 > 5):
                        lam_hb,lam_hr=(temp2[i],temp2[i-1])
        #print('80% wavelength',lam_hb,lam_hr)



        te=np.where((hbeta_flux > 0.9*twe) & (hbeta_flux < 1.1*twe))
        temp=hbeta_range[te]

        if(len(temp) < 2):
                print("asymmetry can't be calculated at 20 %")
                continue

        for p in range(len(temp)):



            if(len(temp)==2):
                    lam_lb,lam_lr=(temp[p],temp[p-1])
            else:
                temp2= temp[p]-temp[p-1]
                if (temp2 > 5):
                    lam_lb,lam_lr=(temp[p],temp[p-1])
        #print('20% wavelength',lam_lb,lam_lr)
        ai=(lam_hb+lam_hr-(lam_lb+lam_lr))/(lam_lr-lam_lb)
        ki= (lam_hr-lam_hb)/ (lam_lr-lam_lb)
        cxh=(lam_hb+lam_hr)/2
        cxl=(lam_lb+lam_lr)/2
        print(ai,ki)
        #print('Asymmetry index:',ai,ki,cxh,cxl)
        #print(q.CalFWHM(0,t),q.CalFWHM(0,s))

        #print(plateid,mjd, fiberid,z,q.SN_ratio_conti,q.CalFWHM(0,t),q.CalFWHM(0,s),ew,area,ratio_n,ratio_b,(integrated/area),q.conti_result[16],blr_size,ai,ki,'\n')
        #fig=plt.figure(figsize=(15,8))
        #plt.plot(hbeta_range,hbeta_flux,color='red',lw=2)
        #plt.plot(q.wave,q.line_flux,color='gray',lw=2)
        #plt.axvline(x=lam_hb,color='blue',lw=1)
        #plt.axvline(x=lam_hr,color='blue',lw=1)
        #plt.axvline(x=lam_lb,color='blue',lw=1)
        #plt.axvline(x=lam_lr,color='blue',lw=1)
        #plt.axvline(x=peakh,color='blue',lw=1)
        #plt.axhline(y=max_flux,color='blue',lw=1)
        #plt.axhline(y=0.75*max_flux,linestyle='--',color='blue',lw=1)
        #plt.axhline(y=0.25*max_flux,linestyle='--',color='blue',lw=1)
        #plt.text(4800,max_flux-20,self.assertNotIn(needle, haystack, 'message'))

        #plt.xlim(4800,4920)
        #plt.ylim(0,max_flux+2)
        #plt.xlabel(r'$\rm Rest \, Wavelength$ ($\rm \AA$)',fontsize = 20)
        #plt.ylabel(r'$\rm f_{\lambda}$ ($\rm 10^{-17} erg\;s^{-1}\;cm^{-2}\;\AA^{-1}$)',fontsize = 20)
        #plt.show()



        '''
        print('H_beta broad component FWHM:',q.CalFWHM(0,t))

    # print('Maximum flux wavelength:',hbeta_range[k],max_flux)
    # print('0.8 flux wavelength',hbeta_range[eg])#,hbeta_flux[eg],eig)
    # print('0.2 flux wavelength',hbeta_range[te])#,hbeta_flux[te],twe)
        #print('----------------------------------------------------------------')



        plt.plot(hbeta_range,hbeta_flux,color='red',lw=2)
        plt.plot(q.wave,q.line_flux,color='gray',lw=3)


        plt.show()
    '''
