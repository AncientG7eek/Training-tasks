from itertools import product
from judge import get_feedback


def generate_possible_codes(n, k):
    """
    Ganerates all the sequences of given length and number of colours

    Args:
        n(int): length of sequence
        k(int): number of colours considered

    Returns:
        list of tuples: all sequences 

    """
    
    return list(product(range(1, k + 1), repeat=n))



def filter_possible_codes(possible_codes, current_guess, feedback):

    """
    Filters the list of possible codes to retain only those that would result in the given feedback if guessed.

    Args:
        possible_codes (list of tuples): A list of all possible codes.
        current_guess (tuple): The current guess being evaluated.
        feedback (tuple): The feedback from the current guess, represented as (correct_place, wrong_place).

    Returns:
        list of tuples: A filtered list of codes that match the feedback for the given guess.
    """
    
    return list(code for code in possible_codes if get_feedback(current_guess, code) == feedback)



def choose_next_guess(all_possible_answers, possible_answers,guesses):

    """
    Chooses the next guess based on minimizing the size of the largest feedback group in possible answers.

    Args:
        all_possible_answers (list of tuples): A list of all possible codes.
        possible_answers (list of tuples): A filtered list of possible codes that are consistent with the feedback.
        guesses (list of tuples): A list of previous guesses.

    Returns:
        tuple: The next guess.

    Raises:
        Prints "scores empty" if no scores are calculated (indicating an issue with the input data).
    """

    scores = {}

    if len(possible_answers) ==1: return possible_answers[0]

    
    for guess in all_possible_answers:
        if guess not in guesses:
            feedback_counts = {}
            
            for possible_answer in possible_answers:
                
                feedback = get_feedback(guess, possible_answer)
                feedback_counts[feedback] = feedback_counts.get(feedback, 0) + 1


            max_group_size = max(feedback_counts.values(), default=0)
            
            scores[guess] = max_group_size
    
    if not scores: print("scores empty")
    
    return min(scores, key=scores.get)