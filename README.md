Support vector regression model for prediction of degron potency from primary sequence.

Composition features are extracted from a list of potential degron sequences (represented as the one-letter amino acid abbreviation) contained in a .csv file by Sequence_to_Composition_reporting.py
Composition features are used by SVM_reporting.r to predict the protein stability index of sequences (range from 1 - 6, lowest stability to highest).
