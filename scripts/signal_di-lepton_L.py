import ROOT
import numpy as np
import math


electron_mass = 0.000511
muon_mass = 0.10566
Z_boson_mass = 91.187
COME = 91 # center of mass energy for e+e- collision
E_res = 1 # energy resolution for histograms

filename = "ntuple.root"

rootfile = ROOT.TFile.Open(filename)

#keys for plots
key1 = 'lepton_pt'
key2 = 'lepton_eta'
key3 = 'lepton_phi'
key4 = 'lepton_theta'
key5 = 'Di-lepton_inv_mass'
key6 = 'p1'
key7 = 'p2'
key8 = 'Di-lepton_pt'
key9 = 'Absolute_p'
key10 = 'Longitudinal_p'
key11 = 'Accoplanarity'
key12 = 'Acolinearity'
key13 = 'Delta_R'

hist_pixel_x = 800
hist_pixel_y = 600

num_hist_x = 2
num_hist_y = 2

num_hist_xx = 2
num_hist_yy = 1

#used two canvases
# LLL- added additional canvas
c1 = ROOT.TCanvas("signal_kinematics1", "signal_kinematics1", hist_pixel_x * num_hist_x, hist_pixel_y * num_hist_y)
c1.Divide(num_hist_x, num_hist_y)
c2 = ROOT.TCanvas("signal_kinematics2", "signal_kinematics2", hist_pixel_x * num_hist_x, hist_pixel_y * num_hist_y)
c2.Divide(num_hist_x, num_hist_y)

c3 = ROOT.TCanvas("signal_kinematics3", "signal_kinematics3", hist_pixel_x * num_hist_x, hist_pixel_y * num_hist_y)
c3.Divide(num_hist_xx, num_hist_yy)

#created lists of kinematic variables so that maximum and minimum values could be obtained for distribution bounds
pts = []
etas = []
phis = []
for entries in rootfile.electron:
	pts.append(entries.electron_pt)
	etas.append(entries.electron_eta)
	phis.append(entries.electron_phi)
for entries in rootfile.muon:
	pts.append(entries.muon_pt)
	etas.append(entries.muon_eta)
	phis.append(entries.muon_phi)

pt_low, pt_high = pts[np.argmin(pts)], pts[np.argmax(pts)]
eta_low, eta_high = etas[np.argmin(etas)], etas[np.argmax(etas)]
phi_low, phi_high = phis[np.argmin(phis)], phis[np.argmax(phis)]

# pt_hist = ROOT.TH1D(key1, 'lepton_pt;pt;Events', 100, int(pt_low) - 1, int(pt_high) + 1)
eta_hist = ROOT.TH1D(key2, 'Lepton eta;#eta;Events', 100, int(eta_low)- .25, int(eta_high) + .25)
phi_hist = ROOT.TH1D(key3, 'Lepton phi;#phi;Events', 100, int(phi_low) - .5, int(phi_high) + .5)
theta_hist = ROOT.TH1D(key4, 'Lepton theta;#theta;Events', 100, 0, 2)
inv_m_hist = ROOT.TH1D(key5, 'Di-lepton Inv Mass;m inv;Events', 100, 0, COME)
# p1_hist = ROOT.TH1D(key6, 'p1;p1;Events', 100, 0, 60)
# p2_hist = ROOT.TH1D(key7, 'p2;p1;Events', 100, 0, 60)
di_lep_pt_hist = ROOT.TH1D(key8, 'Di-lepton pt;pt;Events', 100, 0, int(pt_high) + 1)
di_lep_p_hist = ROOT.TH1D(key9, 'Di-lepton Absolute p;p;Events', 100, 0, 100)
di_lep_pL_hist = ROOT.TH1D(key10, 'Di-lepton Longitudinal p;pL;Events', 100, -50, 50)
acpl_hist = ROOT.TH1D(key11, 'Accoplanarity;Accoplanarity;Events', 100, -3, 10)
#LLL adding more hist, unsure of binnning for delta R for now
acl_hist = ROOT.TH1D(key12, 'Accolinearity;Accolinearity;Events', 100, -3, 10)
delta_R_hist = ROOT.TH1D(key13, 'Delta R;Delta R;Events', 100, 0, 100)


#setting linecolor and title size for each
histograms = [eta_hist, phi_hist, theta_hist, inv_m_hist, di_lep_pt_hist, di_lep_p_hist, di_lep_pL_hist, acpl_hist, acl_hist, delta_R_hist]
for histogram in histograms:
	histogram.SetLineColor(2)
	histogram.SetTitleSize(0.05, 'x')
	histogram.SetTitleSize(0.05, 'y')

def calculate_theta(eta):
	return np.arctan(2 * math.exp(-eta))

def calculate_p(pt, theta):
	return pt / np.sin(theta)

def calculate_pL(pt, theta):
	return pt * np.cos(theta) / np.sin(theta)

iev_prev = 0
num_of_electrons = 0
electron_pt = []
electron_eta = []
electron_phi = []

for entries in rootfile.electron:
	iev = entries.iev
	# move to next particle within an event
	if iev == iev_prev:
		num_of_electrons += 1
		electron_pt.append(entries.electron_pt)
		electron_eta.append(entries.electron_eta)
		electron_phi.append(entries.electron_phi)

	else:
		# reached end of an event, calculate inv mass
		if num_of_electrons >= 2:
			# we select the highest pt leptons
			## find index of max pt electron
			i_max = np.argmax(electron_pt)
			## change max pt to 0 to find index of 2nd max pt  
			electron_pt_max_removed = electron_pt[:]
			electron_pt_max_removed[i_max] = 0
			i_second_max = np.argmax(electron_pt_max_removed)

			# pt_hist.Fill(electron_pt[i_max])
			# pt_hist.Fill(electron_pt[i_second_max])
			eta_hist.Fill(electron_eta[i_max])
			eta_hist.Fill(electron_eta[i_second_max])
			phi_hist.Fill(electron_phi[i_max])
			phi_hist.Fill(electron_phi[i_second_max])
			theta_hist.Fill(calculate_theta(electron_eta[i_max]))
			theta_hist.Fill(calculate_theta(electron_eta[i_second_max]))
			p1 = calculate_p(electron_pt[i_max], calculate_theta(electron_eta[i_max]))
			p2 = calculate_p(electron_pt[i_second_max], calculate_theta(electron_eta[i_second_max]))
			# p1_hist.Fill(calculate_p(p1)
			# p2_hist.Fill(calculate_p(p2)
			di_lep_p_hist.Fill(p1 + p2)
			pL1 = calculate_pL(electron_pt[i_max], calculate_theta(electron_eta[i_max]))
			pL2 = calculate_pL(electron_pt[i_second_max], calculate_theta(electron_eta[i_second_max]))
			di_lep_pL_hist.Fill(pL1 + pL2)
			acpl_hist.Fill(np.pi - (electron_phi[i_max] - electron_phi[i_second_max]))
            #LLL add acolinearity, delta R 
 			#acl_hist.Fill(np.pi - (electron_eta[i_max] - electron_eta[i_second_max]))
 			delta_R_hist.Fill(np.sqrt(np.square(electron_eta[i_max] - electron_eta[i_second_max]) + np.square(electron_phi[i_max] - electron_phi[i_second_max])))

			lep1 = ROOT.TLorentzVector(0,0,0,0)
			lep2 = ROOT.TLorentzVector(0,0,0,0)
			lep1.SetPtEtaPhiM(electron_pt[i_max], electron_eta[i_max],
							electron_phi[i_max], electron_mass)
			lep2.SetPtEtaPhiM(electron_pt[i_second_max], electron_eta[i_second_max], 
							electron_phi[i_second_max], electron_mass)
			pair = lep1 + lep2

			# calculate inv mass and pt from these 2 leptons
			inv_m_hist.Fill(pair.M())
			di_lep_pt_hist.Fill(pair.Perp())

			
		iev_prev = iev
		num_of_electrons = 1
		electron_pt = [entries.electron_pt]
		electron_eta = [entries.electron_eta]
		electron_phi = [entries.electron_phi]


iev_prev = 0
num_of_muons = 0
muon_pt = []
muon_eta = []
muon_phi = []

for entries in rootfile.muon:
	iev = entries.iev
	# move to next particle within an event
	if iev == iev_prev:
		num_of_muons += 1
		muon_pt.append(entries.muon_pt)
		muon_eta.append(entries.muon_eta)
		muon_phi.append(entries.muon_phi)
	else:
		# reached end of an event, calculate inv mass
		if num_of_muons >= 2:
			# we select the highest pt leptons
			## find index of max pt muon
			i_max = np.argmax(muon_pt)
			## change max pt to 0 to find index of 2nd max pt  
			muon_pt_max_removed = muon_pt[:]
			muon_pt_max_removed[i_max] = 0
			i_second_max = np.argmax(muon_pt_max_removed)

			# pt_hist.Fill(muon_pt[i_max])
			# pt_hist.Fill(muon_pt[i_second_max])
			eta_hist.Fill(muon_eta[i_max])
			eta_hist.Fill(muon_eta[i_second_max])
			phi_hist.Fill(muon_phi[i_max])
			phi_hist.Fill(muon_phi[i_second_max])
			theta_hist.Fill(calculate_theta(muon_eta[i_max]))
			theta_hist.Fill(calculate_theta(muon_eta[i_second_max]))
			p1 = calculate_p(muon_pt[i_max], calculate_theta(muon_eta[i_max]))
			p2 = calculate_p(muon_pt[i_second_max], calculate_theta(muon_eta[i_second_max]))
			# p1_hist.Fill(calculate_p(p1)
			# p2_hist.Fill(calculate_p(p2)
			di_lep_p_hist.Fill(p1 + p2)
			pL1 = calculate_pL(muon_pt[i_max], calculate_theta(muon_eta[i_max]))
			pL2 = calculate_pL(muon_pt[i_second_max], calculate_theta(muon_eta[i_second_max]))
			di_lep_pL_hist.Fill(pL1 + pL2)
			acpl_hist.Fill(np.pi - (muon_phi[i_max] - muon_phi[i_second_max]))
            #LLL add acolinearity, delta R 
 			#acl_hist.Fill((np.pi - muon_eta[i_max]) - muon_eta[i_second_max]))
 			delta_R_hist.Fill(np.sqrt(np.square(muon_eta[i_max] - muon_eta[i_second_max]) + np.square(muon_phi[i_max] - muon_phi[i_second_max])))

			lep1 = ROOT.TLorentzVector(0,0,0,0)
			lep2 = ROOT.TLorentzVector(0,0,0,0)
			lep1.SetPtEtaPhiM(muon_pt[i_max], muon_eta[i_max],
							muon_phi[i_max], muon_mass)
			lep2.SetPtEtaPhiM(muon_pt[i_second_max], muon_eta[i_second_max], 
							muon_phi[i_second_max], muon_mass)
			pair = lep1 + lep2

			# calculate inv mass and pt from these 2 lepton
			inv_m_hist.Fill(pair.M())
			di_lep_pt_hist.Fill(pair.Perp())
			
		iev_prev = iev
		num_of_muons = 1
		muon_pt = [entries.muon_pt]
		muon_eta = [entries.muon_eta]
		muon_phi = [entries.muon_phi]



#below is some code found online in an attempt to print two canvases to one pdf to avoid plots being smushed
c1.cd(1)
acpl_hist.Draw()
c1.cd(2)
eta_hist.Draw()
c1.cd(3)
phi_hist.Draw()
c1.cd(4)
theta_hist.Draw()

c2.cd(1)
di_lep_pt_hist.Draw()
c2.cd(2)
di_lep_p_hist.Draw()
c2.cd(3)
di_lep_pL_hist.Draw()
c2.cd(4)
inv_m_hist.Draw()

c3.cd(1)
acl_hist.Draw()
c3.cd(2)
delta_R_hist.Draw()

c1.Print("signal_one.pdf[");
c2.Print("signal_two.pdf");
c3.Print("signal_three.pdf");
