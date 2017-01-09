# MathBotFace

A Python script tweeting about theorems, lemmas and inequalities extracted from wikipedia. (https://twitter.com/McbotfaceBotty)

# Requirements

0. python 3, pip (for arch, `pacman -S python python-pip`)
1. python library wikipedia (`sudo pip install wikipedia`)
2. To tweet, python library tinyurl (to use `tweet.py`)
3. `git clone http://github.com/debdeep777/MathBotFace && cd MathBotFace`

# Usage
* `python fetchWiki.py` to get the titles from wikipedia (which stores it in .DAT files)
* To print a random article on terminal, `python wikiOut.py` or `python wikOut.py inequalities` (with parameter)
* Available parameters for now: `theorems, lemmas, inequalities`
* The default number of characters for the summary is set to 300. Edit the variable `num` in the file to change it.
* To tweet, fill up the Twitter credentials in tweet.py and run it with or without a parameter (e.g. `python tweet.py lemmas`)

