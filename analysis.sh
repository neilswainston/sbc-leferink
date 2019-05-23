pip install --upgrade pip \
	&& pip install -r requirements.txt --upgrade

export PYTHONPATH=$PYTHONPATH:.

python analysis.py out