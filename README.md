* Electric-field spatial TMS-EEG filtering for individualized time-series estimation of primarily targeted neural assemblies 
  
  
![E-field_based_filtering](https://github.com/user-attachments/assets/fbe49df2-44b3-4232-a213-839d321503f4)



* The presented tool enables the creation of source-based spatial filters, using E-field projections to define regions of interest.
The pipeline computes channel weighting to perform spatial filtering. This weighting can be applied for post-processing your data or for real-time filtering in a closed-loop environment. Additionally, it supports the combination of local ROI filters with distant ROI, based on individual DTI cortical projections from the stimulated region.


* General concept:
Concurrent transcranial magnetic stimulation (TMS)  with electroencephalography (EEG) is a popular technique, however, interpreting TMS-EEG results remains challenging due to the no specificity of EEG sources. We introduce a spatial filtering approach where Electric-field simulations are used as spatial filtering boundaries for subject-individualized data curation. The method enables individualized source reconstruction of time series arising from the primarily targeted structures employing a customized minimum norm-type estimator. To validate the present procedure, sensor level analyses was compared to source filtered signals employing TMS evoked potentials. Results suggest that the spatial filtering procedure might be efficient enhancing local hidden neuronal dynamics and reducing components emerging from distant topographies. 




* Recommended pre-reading: 

Niessen E, Bracco M, Mutanen TP, Robertson EM. An analytical approach to identify indirect multisensory cortical activations elicited by TMS? Brain Stimulation: Basic, Translational, and Clinical Research in Neuromodulation 2021;14:376–8.

Hauk O, Stenroos M. A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT. Hum Brain Mapp 2014;35:1642–53.

Hämäläinen MS, Ilmoniemi RJ. Interpreting magnetic fields of the brain: minimum norm estimates. Med Biol Eng Comput 1994;32:35–42.

Taulu S, Simola J. Spatiotemporal signal space separation method for rejecting nearby interference in MEG measurements. Phys Med Biol 2006;51:1759. 


There is no guarantee of any kind.

Xavier COROMINAS,
Martina BRACCO,
Paris, France, 2024.
