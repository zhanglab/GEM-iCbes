import csv

#for row in csv.reader(open('./2019-06-11-updated-deltaG-prediction', mode='rU'), delimiter='\t'):
for row in csv.reader(open('./2019-07-deltaG-predictions', mode='rU'), delimiter='\t'):
	id, dg = row[1], float(row[2])
	dg_kj = dg/1000
	dg_kcal = dg_kj/4.184
	print('{}\t{}\t{}'.format(id, dg_kcal, '4'))
