from __future__ import print_function
import sys
#import trapcam_analysis_old as t
import trapcam_analysis as t

dir = raw_input("Enter the experiment directory you'd like to analyze (e.g. '2017_10_26'): ")
# dir = sys.argv[1]

print ('')
while True:
    analyze_trap_list = []
    letter = raw_input("Enter a trap letter to analyze: ")
    analyze_trap_list.append('trap_'+letter)
    while True:
        letter = raw_input("Enter another trap letter to analyze, or enter 'go' to start batch analysis: ")
        if letter == 'go':
            break
        else:
            analyze_trap_list.append('trap_'+letter)
    print ('')
    print ('you said you want to analyze: ')
    for an_trap in analyze_trap_list:
        print (an_trap)
    user_go_ahead = raw_input("Are those the traps you'd like to analyze? (y/n) ")
    if user_go_ahead == 'y':
        break
    if user_go_ahead == 'n':
        continue
print ('')

calculate_threshold = False
calculate_final = False
thresh_or_final = raw_input("Do you want to analyze just a subset of frames to determine the best in-trap/on-trap threshold, or do you want to do the final analysis? (threshold/final) ")
if thresh_or_final =='threshold':
    calculate_threshold = True
if thresh_or_final == 'final':
    calculate_final = True

for trap in analyze_trap_list:
    analyzer = t.TrapcamAnalyzer(dir, trap, calculate_threshold, calculate_final)
    analyzer.run()
