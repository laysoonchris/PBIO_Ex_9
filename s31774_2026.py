import random
import csv

def validate_positive_integer(prompt: str, min_value: int=1, max_value: int=100_000) -> int:
    while True:
        raw_input = input(prompt)
        try:
            value = int(raw_input)
            if min_value <= value <= max_value:
                return value
            print(f"Please enter an integer between {min_value} and {max_value}.")
        except ValueError:
            print(f"Please enter an integer between {min_value} and {max_value}.")

def validate_seq_id(prompt: str) -> str:
    while True:
        seq_id = input(prompt).strip()
        if seq_id and all(c.isalnum() or c in ('-', '_') for c in seq_id):
            return seq_id
        print("Please enter a valid sequence ID (without whitespace and empty slots).")
