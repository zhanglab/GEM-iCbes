import csv

covered_list = []
for row in csv.reader(open('./2019-07-myco-deltaG-predicitons.tsv', mode='rU'), delimiter='\t'):
	covered_list.append(row[0])

transporter_list = []
for row in csv.reader(open('./myco-transporter-parameters.tsv', mode='rU'), delimiter='\t'):
	transporter_list.append(row[0])

for row in csv.reader(open('./myco-reaction.list', mode='rU'), delimiter='\t'):
	if row[0] not in covered_list:
		if row[0] in transporter_list:
			print('transporter\t{}'.format(row[0]))
		elif 'OMP_' in row[0]:
			print('transporter\t{}'.format(row[0]))
		else:
			print('{}'.format(row[0]))
