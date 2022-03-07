"# randomised_optimization_assignment_2"
Steps
git clone https://github.com/martinmathew/randomised_optimization_assignment_2.git

cd randomised_optimization_assignment_2

conda env update --prefix ./env --file requirements.yml --prune
conda activate path to env\env

#Install Java 8
#RUN Algorithims
java -classpath lib\ABAGAIL.jar;lib\commons-csv-1.9.0.jar;. RandomizedOPT

Generate Charts
python plotcharts.py

For Neural Network
python credit_card.py

