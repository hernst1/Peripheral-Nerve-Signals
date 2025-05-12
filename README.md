# Peripheral-Nerve-Signals

**Analysis of Peripheral Nerve Signals**

**1 Overview**
This is an analysis of peripheral nerve signals (PNSs) recorded through cuff electrodes from a rat’s sciatic nerve. When used for recording, cuff electrodes can capture sensory afferents in response to external stimuli. You will investigate the possibility of discriminating the afferents of three different
sensory stimuli from PNSs.
Details about data acquisition and experimental setup can be found in [1]. The provided dataset has three PNS recordings: each of which was acquired while the rat was in an anesthetized resting state with intermittent application of one of the following stimuli:
• VF: mechanical stimulus of the plantar skin using a Von Frey (VF) filament
• Flex: proprioceptive stimulus provoked by passively flexing the toes
• Pinch: nociceptive stimulus provoked by pinching the toes using fine forceps
The stimuli were specifically chosen to target three di↵erent functional groups of afferent nerve fibers: tactile mechanoreceptive, proprioceptive, and nociceptive [1].

_**1.1 Data Description**_
You will be provided with the following in the ”data.mat” file:
• fs: sampling rate
• VF: data for the mechincal stimulation case
– VF.signal: contains the time series for the VF/rest condition
– VF.trigger: contains the labels for data points (=0 for rest, 6= 0 for stimulation)
• Flex: data for the proprioceptive stimulation case
– Flex.signal: contains the time series for the Flex/rest condition
– Flex.trigger: contains the labels for data points (=0 for rest, 6= 0 for stimulation)
• Pinch: data for the nociceptive stimulation case
– Pinch.signal: contains the time series for the Pinch/rest condition
– Pinch.trigger: contains the labels for data points (=0 for rest, 6= 0 for stimulation)

**2 Tasks**
_**2.1 Pre-processing**_
a) 4th order Butterworth Bandpass Filter, 800-2200 Hz, non-causal
b) Power Spectral Density

_**2.2 Feature Extraction**_
a) MAV and VAR features
b) Feature Selection via SNR

_**2.3 Classification**_

• Built a classifier based on each of the two features separately to discriminate Rest versus Stimulus (VF, Flex, & Pinch) and Stimulus versus Stimulus combinations. Found the 10-fold cross validation accuracy of each of the classifiers and constructed their confusion matrices to answer the following:
• Is it possible to classify all the stimuli against Rest? Is it possible to classify the stimuli from one another? Which of the stimuli is easiest to detect? Which of the two features resulted in a higher accuracy?
• Does the split of the data into training and testing sets follow a certain structure or is it random? Do you think such split would lead to a fair assessment of your classifier’s generalization ability? Why?

**References**
[1] S. Raspopovic, J. Carpaneto, E. Udina, X. Navarro, and S. Micera, ”On the identification of sensory
information from mixed nerves by using single-channel cu↵ electrodes,” Journal of NeuroEngineering
and Rehabilitation, vol. 7, no. 1, p. 17, 2010/04/27 2010.
3
