import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pickle as pkl

def gaus(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def single_smear(gen_eng,obs_eng,save_string):
    arr_min,arr_max,n_bins = np.min(gen_eng),np.max(gen_eng),10.
    bin_len   = (arr_max-arr_min)/n_bins
    x_,y_,z_ = [],[],[]
    for i in range(int(n_bins)):
        bin_midpoint = arr_min + (i + .5)*bin_len
        inds         = np.array([(arr_min + i* bin_len < x < arr_min + (i+1)* bin_len) for x in gen_eng])
        x_0          = obs_eng[inds]
        y, x_edges   = np.histogram(x_0, density=True,bins=25)
        x            = ((np.roll(x_edges,-1) + x_edges)/2.0)[0:len(x_edges)-1]
        mean,sigma   = np.mean(x), np.std(x)

        popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
        mean,sigma = popt[1],popt[2]
        x_.append(bin_midpoint)
        y_.append(mean)
        z_.append(sigma)
        plt.plot(x,y,'b+:',label='data')
        plt.plot(x,gaus(x,*popt),'ro:',label='fit')
        plt.legend()

        plt.xlabel('EDep')
        plt.ylabel('NEntries')
        plt.savefig('../output/validation/%s_%.0f.png' % (save_string,bin_midpoint))
        plt.clf()
    plt.plot(x_,y_)
    plt.xlabel('Beam Energy')
    plt.ylabel('Mean EDep')
    plt.savefig('../output/validation/%s_fit.png' % (save_string))
    plt.clf()

    plt.plot(x_,np.array(z_)/np.array(y_))
    plt.xlabel('Beam Energy')
    plt.ylabel('Std. Dev. of EDep/EDep')
    plt.savefig('../output/validation/%s_std.png' % (save_string))
    plt.clf()
    pkl.dump({'energy_corrections':(x_,y_,z_)},open('../output/pkl/%s_fits.pkl'%(save_string),'wb') )


def merged_smear(bin_1,bin_2,obs_eng,save_string):
    assert len(bin_1) == len(bin_2)
    bin1_min,bin1_max,n_bin1 = np.min(bin_1),np.max(bin_1),5.
    bin1_len   = (bin1_max-bin1_min)/n_bin1

    bin2_min,bin2_max,n_bin2 = np.min(bin_2),np.max(bin_2),5.
    bin2_len   = (bin2_max-bin2_min)/n_bin2
    x_,y_,u_,v_ = [],[],[],[]
    for i in range(int(n_bin1)):
        u_temp,y_temp = [],[]
        for j in range(int(n_bin2)):
            bin1_mp = bin1_min + (i + .5)* bin1_len
            bin2_mp = bin2_min   + (j + .5) *bin2_len

            inds = np.array([(bin1_min + i* bin1_len < bin_1[x] < bin1_min + (i+1)* bin1_len)*(bin2_min + j* bin2_len < bin_2[x] < bin2_min + (j+1)* bin2_len)
                                     for x in range(len(bin_1))])

            x_0          = obs_eng[inds]

            y, x_edges   = np.histogram(x_0, density=True,bins=25)
            x            = ((np.roll(x_edges,-1) + x_edges)/2.0)[0:len(x_edges)-1]
            mean,sigma   = np.mean(x), np.std(x)
            popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
            mean,sigma = popt[1],popt[2]

            x_.append(bin1_mp)
            y_.append(bin2_mp);y_temp.append(bin2_mp)
            u_.append(mean);u_temp.append(mean)

            v_.append(sigma)

            plt.plot(x,y,'b+:',label='data')
            plt.plot(x,gaus(x,*popt),'ro:',label='fit')
            plt.legend()

            plt.xlabel('EDep')
            plt.ylabel('NEntries')
            plt.savefig('../output/validation/%s_%.0f_%.0f.png' % (save_string,bin1_mp,bin2_mp))
            plt.clf()
        plt.plot(y_temp,u_temp)
        plt.xlabel('Bin 2')
        plt.ylabel('Mean EDep')
        plt.savefig('../output/validation/%s_fit_%.0f.png' % (save_string,bin1_mp))
        plt.clf()

    plt.plot(x_,u_)
    plt.xlabel('Beam Energy')
    plt.ylabel('Mean EDep')
    plt.savefig('../output/validation/%s_fit.png' % (save_string))
    plt.clf()

    plt.plot(x_,np.array(v_)/np.array(u_))
    plt.xlabel('Beam Energy')
    plt.ylabel('Std. Dev. of EDep/EDep')
    plt.savefig('../output/validation/%s_std.png' % (save_string))
    plt.clf()
    pkl.dump({'energy_corrections':(x_,y_,u_,v_)},open('../output/pkl/%s_fits.pkl'%(save_string),'wb') )
