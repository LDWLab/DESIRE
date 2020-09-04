# DESIRE

## Installation

1. Initial set up and prerequisites:

	a. Linux distribution (Native Linux system or <a href="https://docs.microsoft.com/en-us/windows/wsl/install-win10#update-to-wsl-2">WSL</a>)

	After it has been set up make sure it is updated:

	>sudo apt-get update

	>sudo apt-get upgrade
	
	b. <a href=" https://code.visualstudio.com/docs/setup/setup-overview">Visual Studio Code</a> (VS Code)

	c. Follow the first two sections of <a href="https://code.visualstudio.com/docs/setup/setup-overview">this guide</a> to install a virtual python environment for the Django project (at least python3.6).

	d. Install Mafft and python3-dev:

	>sudo apt-get install -y mafft
	
	>sudo apt-get install python3-dev
	
	e. Instal git and get an account on github.com. Contact the current admin to receive access to the project repository,
	as well as user account and password for the MySQL database.

	f. In the Linux startup file (.bashrc or .bash_profile) add these lines:

	>export DJANGO_SECRET_KEY='' (provided by admin)
	
	>export DJANGO_USERNAME='' (MySQL username provided by admin)
	
	>export DJANGO_PASSWORD='' (MySQL password provided by admin)

	g. Set up GaTech VPN with the Cisco AnyConnect Secure Mobility Client. Follow <a href="https://faq.oit.gatech.edu/content/how-do-i-get-started-campus-vpn">these</a> instructions.

2. Clone the <a href="https://github.com/LDWLab/DESIRE.git">project repository</a> in a new folder. Get on the latest development branch.

3. Follow the first two sections of <a href="https://code.visualstudio.com/docs/python/tutorial-django">this guide</a> to install a virtual environment in the DESIRE folder.

	a. Prerequisites
	
	b. Create a project environment for the Django tutorial

4. Using the command line from the root directory of the project run the following commands:

	a. Activate the virtual environment

	>source ./env/bin/activate

	b. Install python requirements

	>pip install -r requirements.txt

	c. Install nodejs and its package manager (npm)

	>sudo apt-get install nodejs npm

	d. Install nodejs requirements (listed in package.json)

	>npm install

	e. Initial set up of nodejs scripts

	>./node_modules/.bin/webpack --config webpack.config.js

	f. (Optional) If developing a nodejs script

	>npm run watch

5. Installing PDB topology viewer for development

	a. cd into the pdb-topology-viewer directory and execute the following (might need sudo)

	> npm link

	b. cd into the root of the project directory (one up from the previous location) and run:

	> npm link pdb-topology-viewer

	c. After editing the typescript of the Topology viewer run the following command **from the pdb-topology-directory!**

	> npm run refresh

6. Open a VS Code folder in the root directory of the project.

7. The debugger should now work while connected to the GaTech VPN.
