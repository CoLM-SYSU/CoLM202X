#!/bin/bash

Help()
{
#DISPLAY help
   echo "!----------------------- Usage ----------------------------------------------!"
   echo "Descripion:Summarize Test results"

   echo "!----------------------------------------------------------------------------!"
   echo 'Syntax: ./SummaryTest.bash -n $TestPath/$TestName [-f $TestLists][-i $Varlist]'
   echo "!----------------------------------------------------------------------------!"
   echo options:
   echo "-n The Path and Name of the test working folder"
   echo '-f (Optional) The list of switches of all test cases, if $TestLists '
   echo '   is absent, use $ROOT/run/script/TestLists as the default test list. '
   echo '-i Specify the summary item of the test restuls:' 
   echo '   1)CreateCase;	2)Compile;	3)Submit_Mksrfdata;'
   echo '   4)Submit_Mkinidata;	5)Submit_Case'
   echo '-h display command information'
}

SummaryTest()
{
   ROOT=`awk '/ROOT/ {print $2}' machine.config`
   TestName=$(basename "$1")
   TestPath=$(dirname "$1")

   if [ -z $ROOT ];then
      echo Error in reading ROOT from machine.config
      exit
   fi
   if [ "$3" == "All" ];then
      Varlist="CreateCase Compile Submit_Mksrfdata Submit_Mkinidata Submit_Case"
   fi

TestCaseLists=$2
nfile=`cat $TestCaseLists|wc -l`
for CaseName in `awk '{print $1}' $TestCaseLists`
do
   echo $CaseName
   for Var in $Varlist
   do
      PASS=`cat $TestPath/$TestName/$CaseName/TestStatus|grep $Var|grep PASS`
      if [ ! -z "$PASS" ];then
         echo "	"$PASS
      else
 	 FAIL=`cat $TestPath/$TestName/$CaseName/TestStatus|grep $Var|grep FAIL`
	 if [ ! -z "$FAIL" ];then
    	    echo "	"$FAIL
	 else
            echo "	"$Var PEND
	 fi
      fi
      if [ "$Var" == "Submit_Case" ];then
         if [ -f $TestPath/$TestName/$CaseName/log ];then
            NAN=`grep -i nan $TestPath/$TestName/$CaseName/log`
            if [ -z "$NAN" ];then
               echo "	"NANCheck PASS
            else
               echo "	"NANCheck FAIL
            fi
         else
            echo "	"NANCheck PEND
         fi
         if [ -f $TestPath/$TestName/$CaseName/log ];then
            BALANCE=`grep -i balance $TestPath/$TestName/$CaseName/log`
            if [ -z "$BALANCE" ];then
               echo "	"balance check PASS
            else
               echo "	"balance error FAIL
            fi
         else
            echo "	"balance check PEND
         fi
      fi
   done
done

}


while getopts ":hn:f:i:" options ;
do
    case $options in
      n) TestName="$OPTARG" ;;
      f) TestCaseList="$OPTARG"  ;;
      i) Varlist="$OPTARG" ;;
      h) Help; exit;;
      *) echo "invalid option: $@";exit ;;
    esac
done

if [ -z "${TestName}" ]; then
   echo
   echo 'Error: "-n" is missing, test name is absent'
   echo
   Help
   exit
else
   if [ -z $TestCaseList ];then
      TestCaseList="TestLists"
   fi
fi

if [ -z "$Varlist" ];then
   Varlist="All"
else
   case $Varlist in
      CreateCase) ;;
      Compile);;
      Submit_Mksrfdata);;
      Submit_Mkinidata);;
      Submit_Case);;
      *) echo "invalid Var option: $Varlist"; exit;;
   esac
fi

SummaryTest $TestName $TestCaseList $Varlist
