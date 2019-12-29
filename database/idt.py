#!/home/rdkibler/.conda/envs/domesticator_py36/bin/python
import os
from pathlib import Path
import time
import json
import requests
import base64
import argparse

user_info_file = os.path.expanduser("~/.domesticator/info.json")
token_file = os.path.expanduser("~/.domesticator/token.json")

def vprint(str,verbose=False,**kwargs):
	if verbose: print(str,**kwargs)

def parse_user_args():
	parser=argparse.ArgumentParser(prog='domesticator', description='The coolest codon optimizer on the block')

	parser.add_argument('seqs', nargs='*', help="supply a space separated DNA sequences to query")
	parser.add_argument("--interactive",default=False,action="store_true")
	parser.add_argument('--store_token',default=False,action="store_true")
	parser.add_argument("--new_token",default=False,action="store_true")
	parser.add_argument("--delete_token",default=False,action="store_true")
	parser.add_argument("--reset_user_info",default=False,action="store_true")
	parser.add_argument("--verbose",default=False,action="store_true")

	return parser.parse_args()

def ask_for_user_data(user_info_file):
	user_info = {}
	idt_url = "https://www.idtdna.com/pages/tools/apidoc"
	msg = "The first time you use the IDT " \
		"complexity function requires you to " \
		"input your IDT username, password, " \
		"API client ID and API client secret. "\
		"These will be stored securely in " \
		f"'{os.path.abspath(user_info_file)}', " \
		"which only you have access to. Before " \
		"you begin, go to {idt_url} and follow " \
		"all the directions under the 'Get " \
		"access to the API' header. The client " \
		"ID and client description can be " \
		"anything. The client secret will be " \
		"generated for you."

	print(msg)
	print("1) Please enter your IDT account username: ")
	user_info["username"] = input()
	print("2) Please enter your IDT account password: ")
	user_info["password"] = input()
	print("3) Please enter you API client ID: ")
	user_info["ID"] = input()
	print("4) Please enter your API secret: ")
	user_info["secret"] = input()

	return user_info

def get_user_info(user_info_file):
	if os.path.exists(user_info_file):
		with open(user_info_file,'r') as f:
			user_info = json.load(f)

	else:
		#we have to set things up for the first time
		os.makedirs(os.path.dirname(user_info_file),exist_ok=True)

		user_info = ask_for_user_data(user_info_file)

		#do it this weird way just to make sure nobody can see your private stuff
		#simply creates an empty file (like touch)
		Path(user_info_file).touch()
		#make it so only the user can acces this file
		os.chmod(user_info_file,0o600)
		#now write the secret stuff ;)
		with open(user_info_file,'w+') as f:
			json.dump(user_info,f)

	return user_info

def delete_stored_token(token_file):
	if os.path.exists(token_file):
		os.remove(token_file)

def store_token(token,token_file):
	delete_stored_token(token_file)
	#do it this weird way just to make sure nobody can see your private stuff
	#creates an empty file
	Path(user_info_file).touch()
	#make it so only the user can acces this file
	os.chmod(user_info_file,0o600)
	#now write the secret stuff ;)
	with open(token_file,'w+') as f:
		json.dump(token,f)


def get_new_token(user_info, verbose=False):
	vprint("getting new token", verbose)
	client_info = user_info["ID"]+":"+user_info["secret"]
	byte_encoded_client_info = client_info.encode("utf-8")
	base64_client_info = base64.urlsafe_b64encode(byte_encoded_client_info).decode('utf8')

	url = "https://www.idtdna.com/Identityserver/connect/token"

	payload = f'grant_type=password&username={user_info["username"]}&password={user_info["password"]}&scope=test'
	headers = {
		'Content-Type': 'application/x-www-form-urlencoded',
		'Authorization': f'Basic {base64_client_info}'
	}
	tries = 5
	while(tries > 0):
		response = requests.request("POST", url, headers=headers, data = payload)

		response_dict = json.loads(response.text)
		tries =- 1
		if "access_token" in response_dict:
			vprint("token acquired",verbose)
			break
		else:
			vprint(f"ERROR: {response_dict['Message']}\nTrying again in 5 seconds...",verbose)
			time.sleep(5)
	if "access_token" not in response_dict:
		raise RuntimeError("Could not get access token")
	return response_dict

def get_stored_token(token_file):
	with open(token_file,'r') as f:
		return json.load(f)

def get_token(token_file, user_info,verbose=False):
	get_token_flag = False
	if os.path.exists(token_file):
                vprint(f"using token stored at {token_file}",verbose)
                modified_time = os.path.getmtime(token_file)
                current_time = time.time()
                token = get_stored_token(token_file)
                if current_time - modified_time > token["expires_in"]:
                        vprint("token expired",verbose)
			get_token_flag = True
        else:
                vprint(f"no file found at {token_file}",verbose)
		get_token_flag = True

	if get_token_flag:
                token = get_new_token(user_info)
                vprint(f"storing token at {token_file}",verbose)
                store_token(token,token_file)
	return token

def query_complexity(seq, token,verbose=False):
	vprint("querying complexity",verbose)
	url = "https://www.idtdna.com/api/complexities/screengBlockSequences"

	payload = f'[{{"Name":"My gBlock","Sequence":"{seq}"}}]'
	headers = {
		'Content-Type': 'application/json',
		'Authorization': f'Bearer {token}'
	}

	response = requests.request("POST", url, headers=headers, data = payload)

	return json.loads(response.text)

def main():
	args = parse_user_args()

	if args.reset_user_info:
		vprint(f"deleting {user_info_file}",args.verbose)
		os.remove(user_info_file)

	vprint("getting user info",args.verbose)
	user_info = get_user_info(user_info_file)

	if args.delete_token:
		vprint(f"deleting token stored at {token_file}",args.verbose)
		delete_stored_token(token_file)
		exit("token deleted")

	if args.new_token:
		token = get_new_token(user_info, args.verbose)
	else:
		token = get_token(token_file, user_info, args.verbose)

	if args.store_token:
		vprint(f"storing token at {token_file}",args.verbose)
		store_token(token,token_file)

	if args.interactive:
		while(True):
			print("input a DNA sequence or hit return without typing anything else to exit")
			seq = input()
			if seq == "":
				break
			response = query_complexity(seq,token["access_token"])
			if len(response[0]) == 0:
				print("no issue!")
			score_sum = 0
			for issue in response[0]:
				vprint(issue,args.verbose)
				print(issue["Score"],issue["Name"])
				score_sum += issue["Score"]
			print(f"Total Score: {score_sum}")


	else:
		for seq in args.seqs:
			response = query_complexity(seq,token["access_token"])
			if len(response[0]) == 0:
				print("no issue!")
			score_sum = 0
			for issue in response[0]:
				vprint(issue,args.verbose)
				print(issue["Score"],issue["Name"])
				score_sum += issue["Score"]
			print(f"Total Score: {score_sum}")



if __name__ == "__main__":
	main()
"""
Accepted - Moderate Complexity (Scores between 7 and 19)

Some complexities exist that may interfere with or delay manufacturing. If it is possible to reduce these complexities please do so, otherwise we will attempt this order.
"""
