#!/bin/bash
 

time0=$(date +%s%N | cut -b1-13) #start time
 
#deletes first 2 characters of a line:
cut_it() {
    echo $1 | cut -c 3-
}
 
echo "=== delete old files"
rm -rf src
 
echo "=== COPY SRC "
cp -r ../../src .
cp mainpage.txt src/.
 
echo "============================================================================================================================" > missing_comment_lines.txt
 
echo "=== delete old version of documentation"
rm -rf doxygen

echo "=== Delete empty comments"
FILELIST=`find . -name *.f90`
 for file in $FILELIST
do
  perl -pi -e 's/! \?/!  /g' $file
done 


echo "=== delete INTERFACE/END INTERFACE lines and MODULE PROCEDURE lines"
for file in $FILELIST
do
	perl -pi -e 's/END INTERFACE/!END INTERFACE/g' $file
	perl -pi -e 's/INTERFACE/!INTERFACE/g' $file
	perl -pi -e 's/MODULE PROCEDURE/!MODULE PROCEDURE/g' $file
done
 
echo "=== Comments for input and output variables"
for file in $FILELIST
do
  perl -pi -e 's/\(IN\)(.*):: (.*?)! (.*)/\(IN\)$1:: $2! $3\n!> \\param[in] $2\n!> \\verbatim\n!>$3\n!> \\endverbatim/g' $file
  perl -pi -e 's/\(OUT\)(.*):: (.*?)! (.*)/\(OUT\)$1:: $2! $3\n!> \\param[out] $2\n!> \\verbatim\n!>$3\n!> \\endverbatim/g' $file
  perl -pi -e 's/\(INOUT\)(.*):: (.*?)! (.*)/\(INOUT\)$1:: $2! $3\n!> \\param[inout] $2\n!> \\verbatim\n!>$3\n!> \\endverbatim/g' $file
done

echo "=== Comments for MODULE, SUBROUTINE, RECURSIVE SUBROUTINE, FUNCTION, PURE FUNCTION and PROGRAM"
 
name=(SUBROUTINE RECURSIVE MODULE FUNCTION PURE PROGRAM) 
name2=('' ' SUBROUTINE' '' '' ' FUNCTION' '') 
Len_name=${#name[@]} #Number of names
 
#for time calculating:
time_all=0
y=1
z=0
for (( k=0; k<${Len_name}; k++ )); 
do 
  FILELIST2=($(grep -l "^"${name[k]} $FILELIST | sed 's/:.*//')) 
  z=$[$z+${#FILELIST2[@]}]
done
 
for (( k=0; k<${Len_name}; k++ )); 
do 
   
  echo "=== Comments for "${name[k]}${name2[k]}
  echo "====================================================" 
   
  FILELIST2=($(grep -l "^"${name[k]} $FILELIST | sed 's/:.*//')) 
  Length=${#FILELIST2[@]} #Number of matching files
   
  for (( j=0; j<${Length}; j++ ));
  do
   
    time1=$(date +%s%N | cut -b1-13) #time at begin of loop
     
    echo "=== File "$[j+1]" of "$Length" for "${name[k]}${name2[k]}" comments" 
    echo "=== Add empty Doxygen comments in "${FILELIST2[j]}" for "${name[k]}${name2[k]}" comments" 
    perl -pi -e 's/^'${name[k]}'(.*)\n(.*)/!> \\par Purpose:\n!> \\verbatim\n!> \n!> \\endverbatim\n'${name[k]}'$1\n/g' ${FILELIST2[j]}
     
    echo "=== Read "${FILELIST2[j]}" line by line"
    line=()
    typeset -i i=0
    while read line[$i]
    do
      i=i+1
    done < ${FILELIST2[j]}
     
    Len=${#line[@]} #Number of lines
     
    echo "=== Copy "${name[k]}${name2[k]}" comments to doxygen comment lines in "${FILELIST2[j]}
     
    for (( i=${Len}; i>=0; i-- ));
    do 
      if echo ${line[i]} | grep "^"${name[k]} >/dev/null
      then 
      
      	and=0
      	if echo ${line[i]} | grep "&$" >/dev/null;then #true if variables are in two lines
      	  and=1
      	fi
      
        if echo ${line[i+1+and]} | grep "!============================================================================================================================" >/dev/null ;then
          if echo ${line[i+2+and]} | grep "!============================================================================================================================" >/dev/null ;then
            echo "Missing comment in line: "$[i+1]". File: "${FILELIST2[j]}"." |cat - missing_comment_lines.txt > /tmp/out && mv /tmp/out missing_comment_lines.txt
          else
            for (( l=3; l<99; l++ ));
            do  
              if echo ${line[i+l]} | grep "!============================================================================================================================" >/dev/null;then
                if echo ${line[i+l-1]} | grep "^!" >/dev/null; then
                  x=0
                  for (( m=0; m< (${l}-${and}-3); m++ )); do  
                    sed -i -e $[i-1]c"!>\n!>" ${FILELIST2[j]} #adds more empty lines if comment is not a one-liner
                    x=$[$x+1]
                  done
                  for (( n=0; n<=(l-and-3); n++ )); do  
                    sed -i -e  $[i-1+n]c"!> $(cut_it "${line[i+l+n-1-x]}")" ${FILELIST2[j]} #copies comment to doxygen comment line(s)
                  done
                else
                  echo "Missing comment in line: "$[i+l-1]". File: "${FILELIST2[j]}"." |cat - missing_comment_lines.txt > /tmp/out && mv /tmp/out missing_comment_lines.txt
                fi
                break
              else
                continue
              fi
            done  
          fi
        fi
      fi
    done
     
    #for time calculating:
    time2=$(date +%s%N | cut -b1-13) #time at end of loop
    time_d=$(($time2 - $time1)) 
    time_all=$(($time_all+$time_d)) 
    time_e=$(($time_all/$y)) 
    echo "=== Remaining time:" $((($time_e * ($z - $y))/1000)) "seconds" #time
    echo "=== Progress: " $((100*$y / $z))"%" #time
    y=$[$y+1] #time
     
    echo "====================================================" 
   
  done
 
done
 
echo "=== create documentation with modified source files"
doxygen doxyconfig

echo "=== edit files.html and src-html files"
cd doxygen/html/

k=1

cp files.html files.html.bak

line=()
typeset -i i=0
while read line[$i]
do
  i=i+1
done < files.html

Len=${#line[@]} #Number of lines

for (( i=${Len}; i>=0; i-- ));
do 
  if echo ${line[i]} | grep "_8f90.html" >/dev/null
   		then echo $k >/dev/null
  		k=$[$k+1]
  		perl -pi -e 's/(.*)class="el" href="(.*)_8f90.html"(.*)/$2_8f90.html/g' files.html 
  fi
done

line2=()
typeset -i j=0
while read line2[$j]
do
  j=j+1
done < files.html

Len2=${#line2[@]} #Number of lines

for (( j=${Len2}; j>=0; j-- ));
do 
if echo ${line2[j]} | grep "_8f90.html" >/dev/null
 	then echo ${line2[j]} >/dev/null
	
	linex=()
	typeset -i x=0
	while read linex[$x]
	do
 		 x=x+1
	done < ${line2[j]}

	perl -pi -e 's/(.*)classmod__(.*).html"(.*)/classmod__$2.html/g' ${line2[j]}

	liney=()
	typeset -i y=0
	while read liney[$y]
	do
 		 y=y+1
	done < ${line2[j]}
	Leny=${#liney[@]} #Number of lines

	for (( y=${Leny}; y>=0; y-- ));
	do 
	if echo ${liney[y]} | grep "classmod__" >/dev/null
 		then echo ${liney[y]} >/dev/null
		echo "COPY "${liney[y]} " TO " ${line2[j]}
		cp ${liney[y]} ${line2[j]}
	fi
	done
	
fi
done

mv files.html.bak files.html

perl -pi -e 's/!END !INTERFACE/END INTERFACE/g' *_8f90_source.html
perl -pi -e 's/!INTERFACE/INTERFACE/g' *_8f90_source.html
perl -pi -e 's/!MODULE PROCEDURE/MODULE PROCEDURE/g' *_8f90_source.html

cd ../..

time3=$(date +%s%N | cut -b1-13) #finish time 
 
echo "=== DONE. TOTAL TIME: "$((($time3 - $time0)/1000))" seconds"

cp HOPR*.png doxygen/html/.
 
line=()
typeset -i i=0
while read line[$i]
do
  i=i+1
done < missing_comment_lines.txt
if echo ${line[0]} | grep "============================================================================================================================" >/dev/null;then
  rm missing_comment_lines.txt
else 
  echo "=== Script found missing comment lines!!! Results in HOPR/doxygen/missing_comment_lines.txt"
fi
