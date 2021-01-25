import ROOT
import numpy as np


electron_mass = 0.000511
muon_mass = 0.10566
COME = 91 # center of mass energy for e+e- collision
E_res = 1 # energy resolution for histograms

filename1 = "ntuple.root"
filename2 = "ntuple_bkg.root"

rootfile1 = ROOT.TFile.Open(filename1)
rootfile2 = ROOT.TFile.Open(filename2)
key1 = 'pt_s'
key2 = 'pt_2f_bkg'
key3 = 'inv_m_s'
key4 = 'inv_m_2f_bkg'

hist_pixel_x = 800
hist_pixel_y = 600
num_hist_x = 2
num_hist_y = 1


canvas = ROOT.TCanvas("pt", "pt", hist_pixel_x * num_hist_x, hist_pixel_y * num_hist_y)
canvas.Divide(num_hist_x, num_hist_y)

pts_s = [entries.electron_pt for entries in rootfile1.electron]
pts_s += [entries.muon_pt for entries in rootfile1.muon]
pts_b = [entries.electron_pt for entries in rootfile2.electron]
pts_b += [entries.muon_pt for entries in rootfile2.muon]

pt_s_low, pt_s_high = pts_s[np.argmin(pts_s)], pts_s[np.argmax(pts_s)]
pt_b_low, pt_b_high = pts_b[np.argmin(pts_b)], pts_b[np.argmax(pts_b)]

pt_s = ROOT.TH1D(key1, 'Di-lepton pt;pt;Events/GeV', 100, 0, max(pt_s_high, pt_b_high))
pt_b = ROOT.TH1D(key2, 'Di-lepton pt;pt;Events/GeV', 100, 0, max(pt_s_high, pt_b_high))
inv_m_s = ROOT.TH1D(key3, 'Di-lepton Inv Mass;m inv;Events/GeV', 100, 30, COME)
inv_m_b = ROOT.TH1D(key4, 'Di-lepton Inv Mass;m inv;Events/GeV', 100, 30, COME)

sig_events = 0
bkg_events = 0

for file in [rootfile1, rootfile2]:

	iev_prev = 0
	num_of_electrons = 0
	electron_pt = []
	electron_eta = []
	electron_phi = []

	for entries in file.electron:
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
				lep1 = ROOT.TLorentzVector(0,0,0,0)
				lep2 = ROOT.TLorentzVector(0,0,0,0)
				lep1.SetPtEtaPhiM(electron_pt[i_max], electron_eta[i_max],
								electron_phi[i_max], electron_mass)
				lep2.SetPtEtaPhiM(electron_pt[i_second_max], electron_eta[i_second_max], 
								electron_phi[i_second_max], electron_mass)
				pair = lep1 + lep2
				inv_m = pair.M()
				pt = pair.Pt()
				if inv_m > 30:
					if file == rootfile1:
						pt_s.Fill(pt)
						inv_m_s.Fill(inv_m)
						sig_events += 1
					else:
						pt_b.Fill(pt)
						inv_m_b.Fill(inv_m)
						bkg_events += 1

				# calculate inv mass from these 2 leptons
				# inv_m.Fill((lep1+lep2).M())
				# recoil_m.Fill(COME ** 2 + ((lep1+lep2).M()) ** 2 - 2 * COME * (lep1+lep2).E())
			iev_prev = iev
			num_of_electrons = 1
			electron_pt = [entries.electron_pt]
			electron_eta = [entries.electron_eta]
			electron_phi = [entries.electron_phi]

for file in [rootfile1, rootfile2]:

	iev_prev = 0
	num_of_muons = 0
	muon_pt = []
	muon_eta = []
	muon_phi = []

	for entries in file.muon:
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
				lep1 = ROOT.TLorentzVector(0,0,0,0)
				lep2 = ROOT.TLorentzVector(0,0,0,0)
				lep1.SetPtEtaPhiM(muon_pt[i_max], muon_eta[i_max],
								muon_phi[i_max], muon_mass)
				lep2.SetPtEtaPhiM(muon_pt[i_second_max], muon_eta[i_second_max], 
								muon_phi[i_second_max], muon_mass)
				pair = lep1 + lep2
				inv_m = pair.M()
				pt = pair.Pt()
				if inv_m > 30:
					if file == rootfile1:
						pt_s.Fill(pt)
						inv_m_s.Fill(inv_m)
						sig_events += 1
					else:
						pt_b.Fill(pt)
						inv_m_b.Fill(inv_m)
						bkg_events += 1

				# calculate inv mass from these 2 lepton
				# inv_m.Fill((lep1+lep2).M())
				# recoil_m.Fill(COME ** 2 + ((lep1+lep2).M()) ** 2 - 2 * COME * (lep1+lep2).E())
			iev_prev = iev
			num_of_muons = 1
			muon_pt = [entries.muon_pt]
			muon_eta = [entries.muon_eta]
			muon_phi = [entries.muon_phi]


pt_s.Scale(1/pt_s.Integral())
pt_b.Scale(1/pt_b.Integral())
inv_m_s.Scale(1/inv_m_s.Integral())
inv_m_b.Scale(1/inv_m_b.Integral())

pt_s.SetLineColor(ROOT.kRed)
pt_s.SetMinimum(0)
pt_b.SetLineColor(ROOT.kBlue)
pt_b.SetMinimum(0)

inv_m_s.SetLineColor(ROOT.kRed)
inv_m_s.SetMinimum(0)
inv_m_b.SetLineColor(ROOT.kBlue)
inv_m_b.SetMinimum(0)

canvas.cd(1)
inv_m_s.Draw("HIST")
inv_m_b.Draw("HIST SAME")

canvas.cd(2)
pt_s.Draw("HIST")
pt_b.Draw("HIST SAME")


canvas.SaveAs("inv_mass_pt_cut_2fbkg.pdf")
