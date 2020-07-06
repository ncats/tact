#!/usr/bin/python

# Ref: https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
# Ref: https://stackabuse.com/read-a-file-line-by-line-in-python/


data = open ('../input/Malaria_TAUC_sliding_window_concat.tab', 'r')

line = data.readline ()
newline = ''
cnt = 0

response = ''
cell_line = ''
header = 'cell.line'

while line:
	#print (line.strip())


	if cnt == 0:
		line = line.replace (' ', '_')
		tmp = line.strip().split('\t')
		tmp = tmp[0:]
		for i in range(1, 12):
			for j in range(1, 13 - i):
				for h in tmp:
					header += '\t' + h + '_' + str(i) + '_' + str(j)
		print (header)

	if cnt > 0:


		tmp = line.strip().split('\t')

		cell_line = tmp[0].replace('_1_1', '').replace('_', ' ').replace('KH', '163-KH').replace('KH1-001', 'KH1-001-RME').replace('KH1-060', 'KH1-060-RME')
		#cell_line = tmp[0].split('_')[0]


		response = tmp [1:]
		
		if cnt == 1:
			newline = cell_line
		
		for r in response:
			newline += '\t' + r


	cnt += 1
	line = data.readline ()
	
	if cnt == 67:
		print (newline)

		cnt = 1
		newline = ''

data.close()
