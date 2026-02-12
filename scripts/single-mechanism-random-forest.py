import pandas
import numpy as np

import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

pathwaysDataframe = pandas.read_csv('Morgan Difference Pathways Dataset.csv')

columns2drop = ['Reaction Mechanism', 'Mechanism Class', 'Mechanism Bond Changes']
morganDiffVectors = pathwaysDataframe.drop(columns2drop, axis=1).values # morgan difference fingerprint features

mechanismClasses = pathwaysDataframe['Mechanism Class'].values # mechanism class labels

morganDiffVectorsTrain, morganDiffVectorsTest, mechanismClassesTrain, mechanismClassesTest = train_test_split(
    morganDiffVectors, mechanismClasses, test_size=0.2, random_state=42, stratify=mechanismClasses
)

mechanismClassifier = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)

mechanismClassifier.fit(morganDiffVectorsTrain, mechanismClassesTrain)

mechanismClassesPredicter = mechanismClassifier.predict(morganDiffVectorsTest)
mechanismClassificationAccuracy = accuracy_score(mechanismClassesTest, mechanismClassesPredicter)

print("Mechanism Classification Accuracy:", mechanismClassificationAccuracy)


