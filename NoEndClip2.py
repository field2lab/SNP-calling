# NoEndClip SCRIPT
# Recover soft clipped bases at ends
#edit CIGAR, PNEXT, and TLEN variables

# Module Imports
import sys, os, time
import re
import argparse
start_time = time.time()

parser = argparse.ArgumentParser(description='Edit soft clipped base')
parser.add_argument('-isam', '--input_sam', required=True, type=str)
parser.add_argument('-osam', '--output_sam', required=True, type=str)
#set_inference_options(parser)
args = parser.parse_args()

path_to_old_sam = args.input_sam
path_to_new_sam = args.output_sam

# ----------------------------------------------------------------------------------------------------
# Recover SOFTCLIP line if any
def Modify_fields(InputRow):
	leftclipNum = 0
	rightclipNum = 0
	Outrow = [InputRow,leftclipNum,rightclipNum]
	#if any('ZA:Z' in k for k in InputList):
	#if all(x in InputRow for x in ['ZA:Z','1S']):
	if '1S' in InputRow:
		InputList = InputRow.split(sep = '\t')
		if (('ZA:Z' not in InputRow) or (len([x for x in InputList if 'ZA:Z' in x][0]) < 11)):


			#----------edit R1:
				# soft clip left end on R1
				if (InputList[5][:2] == '1S' and int(InputList[8]) > 0 and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[8] = str(int(InputList[8])+1)
					leftclipNum += 1
				# soft clip right end on R1
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) > 0):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					rightclipNum += 1
			#----------edit R2:
				# soft clip left end on R2
				if (InputList[5][:2] == '1S' and int(InputList[8]) < 0  and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					leftclipNum += 1
				# soft clip right end on R2
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) < 0):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[8] = str(int(InputList[8])-1)
					rightclipNum += 1
			# soft clip left end on R1/R2 and pair strand unmapped
				if (InputList[5][:2] == '1S' and int(InputList[8]) == 0 and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[7] = str(int(InputList[7])-1)
					leftclipNum += 1
			# soft clip right end on R1/R2 and pair strand unmapped
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) == 0):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					rightclipNum += 1

		elif len([x for x in InputList if 'ZA:Z' in x][0]) >= 11:
				ZAstr = [x for x in InputList if 'ZA:Z' in x][0]
				ZAcigar = ZAstr.split(';')[11]
				if ZAstr.split(';')[0] == "ZA:Z:<&":
					ZAcigar = ZAstr.split(';')[5]
				#print(ZAcigar)
			#----------edit R1:
				# soft clip left end on R1 and no clip on R2
				if (InputList[5][:2] == '1S' and int(InputList[8]) > 0 and ("S" not in ZAcigar) and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[8] = str(int(InputList[8])+1)
					leftclipNum += 1
				# soft clip right end on R1 and no clip on R2
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) > 0 and ("S" not in ZAcigar)):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					rightclipNum += 1
				# soft clip left end on R1 and R2 has left soft clip
				if (InputList[5][:2] == '1S' and int(InputList[8]) > 0 and ZAcigar[:2] == '1S' and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[7] = str(int(InputList[7])-1)
					InputList[8] = str(int(InputList[8])+1)
					leftclipNum += 1
				# soft clip left end on R1 and R2 has right soft clip
				if (InputList[5][:2] == '1S' and int(InputList[8]) > 0 and ZAcigar[len(InputList[5])-2:] == '1S' and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[8] = str(int(InputList[8])+2)
					leftclipNum += 1
				# soft clip right end on R1 and R2 has left soft clip
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) > 0 and ZAcigar[:2] == '1S'):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[7] = str(int(InputList[7])-1)
					rightclipNum += 1
				# soft clip right end on R1 and R2 has right soft clip
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) > 0 and ZAcigar[len(ZAcigar)-2:] == '1S'):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[8] = str(int(InputList[8])+1)
					rightclipNum += 1
				# no clip on R1 and R2 has left soft clip
				if (("S" not in InputList[5]) and int(InputList[8]) > 0 and ZAcigar[:2] == '1S'):
					InputList[7] = str(int(InputList[7])-1)
				# no clip on R1 and R2 has right soft clip
				if (("S" not in InputList[5]) and int(InputList[8]) > 0 and ZAcigar[len(ZAcigar)-2:] == '1S'):
					InputList[8] = str(int(InputList[8])+1)

			#----------edit R2:
				# soft clip left end on R2 and no clip on R1
				if (InputList[5][:2] == '1S' and int(InputList[8]) < 0 and 'S' not in ZAcigar):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					leftclipNum += 1
				# soft clip right end on R2 and no clip on R1
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) < 0 and 'S' not in ZAcigar):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[8] = str(int(InputList[8])-1)
					rightclipNum += 1
				# soft clip left end on R2 and R1 has left soft clip
				if (InputList[5][:2] == '1S' and int(InputList[8]) < 0 and ZAcigar[:2] == '1S' and int(InputList[7]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[7] = str(int(InputList[7])-1)
					InputList[8] = str(int(InputList[8])-1)
					leftclipNum += 1
				# soft clip left end on R2 and R1 has right soft clip
				if (InputList[5][:2] == '1S' and int(InputList[8]) < 0 and ZAcigar[len(ZAcigar)-2:] == '1S'):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					leftclipNum += 1
				# soft clip right end on R2 and R1 has left soft clip
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) < 0 and ZAcigar[:2] == '1S'):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[7] = str(int(InputList[7])-1)
					InputList[8] = str(int(InputList[8])-2)
					rightclipNum += 1
				# soft clip right end on R2 and R1 has right soft clip
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) < 0 and ZAcigar[len(ZAcigar)-3:] == 'M1S'):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					InputList[8] = str(int(InputList[8])-1)
					rightclipNum += 1
				# no clip on R2 and R1 has left soft clip
				if (("S" not in InputList[5]) and int(InputList[8]) < 0 and ZAcigar[:2] == '1S' and int(InputList[7]) != 1):
					InputList[7] = str(int(InputList[7])-1)
					InputList[8] = str(int(InputList[8])-1)

			# soft clip left end on R1/R2 and pair strand unmapped
				if (InputList[5][:2] == '1S' and int(InputList[8]) == 0 and int(InputList[3]) != 1):
					matchbase = ''.join(re.findall(r"1S(.*?)M",InputList[5]))
					#print(matchbase)
					leftstr = str(int(matchbase)+1)+'M'
					#print(leftstr)
					InputList[5] = leftstr + InputList[5].replace(("1S"+matchbase+'M'), "")
					InputList[3] = str(int(InputList[3])-1)
					InputList[7] = str(int(InputList[7])-1)
					leftclipNum += 1
			# soft clip right end on R1/R2 and pair strand unmapped
				if (InputList[5][len(InputList[5])-3:] == 'M1S' and int(InputList[8]) == 0):
					matchbase = ''.join(re.findall(r"(\d+)M1S",InputList[5]))
					#print(matchbase)
					rightstr = str(int(matchbase)+1)+'M'
					#print(rightstr)
					InputList[5] = InputList[5].replace((matchbase+'M1S'), "")+rightstr
					rightclipNum += 1
			# unmapped strand but paired R1/R2 has soft clip left
				if (InputList[5] == '*' and int(InputList[3]) != 1 and ZAcigar[:2] == '1S'):
					InputList[3] = str(int(InputList[3])-1)
					InputList[7] = str(int(InputList[7])-1)

			# Return Output
		Outrow = ['	'.join(InputList),leftclipNum,rightclipNum]
	return Outrow

# ----------------------------------------------------------------------------------------------------
with open(path_to_old_sam, "r") as old_sam:
	with open(path_to_new_sam, "w") as new_sam:
		leftclipTotal = 0
		rightclipTotal = 0
		for i,line in enumerate(old_sam):
			lineparse = Modify_fields(line)
			leftclipTotal += lineparse[1]
			rightclipTotal += lineparse[2]
			new_sam.writelines(lineparse[0])
	new_sam.close()
old_sam.close()
print("total leftclips recovered: " + str(leftclipTotal))
print("total rightclips recovered: " + str(rightclipTotal))
print("--- %s seconds ---" % round((time.time() - start_time),4))
