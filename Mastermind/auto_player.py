import Interactive_Module as IM
import Auto_Module as AM
from judge import get_feedback



hidden_code = IM.query_seq("\nEnter a sequence you want the computer to guess:\n>", "Enter a sequence of space separated numbers:\n")
n = len(hidden_code)
k = IM.query_num("\nEnter the number of colours considered:\n>", "Enter a number:\n")

print("\nComputer is about to take up the challange!\n")

possible_answers=None
guesses=None
all_codes = AM.generate_possible_codes(n, k)


while True:
    if possible_answers is None:
        possible_answers = all_codes
    if guesses is None:
        guesses = []


    if not guesses:
        current_guess = tuple(1 for _ in range(n // 2)) + tuple(
            2 for _ in range(n - n // 2)
        )
        guesses.append(current_guess)
    else:
        current_guess = guesses[-1]


    feedback = get_feedback(current_guess, hidden_code)

    print(f"Round {len(guesses)}. Guess: {current_guess}, Feedback: {feedback}, Number of possible answers left: {len(possible_answers)}")

    # if the code is found
    if feedback == (n, 0):
        print(f"Code found: {current_guess}")
        print(f"Guesses: {guesses}")
        print(f"Done! Computer broke your code in {len(guesses)} attempts.")
        break

    # narrowing the possible answers space based on new feedback
    possible_answers = AM.filter_possible_codes(possible_answers, current_guess, feedback)
    

    next_guess = AM.choose_next_guess(
        all_codes, possible_answers, guesses)

    guesses.append(next_guess)
