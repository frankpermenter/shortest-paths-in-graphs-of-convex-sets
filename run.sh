set -e
export PYTHONPATH="${PYTHONPATH}:$PWD/spp"
#export PYTHONPATH="${PYTHONPATH}:/opt/drake/lib/python3.6/site-packages/"
export PYTHONPATH="${PYTHONPATH}:/home/frank/drake-build/install/lib/python3.6/site-packages/"

#set path=$PWD
#cd $path
python3 run.py
#python3 simple.py

