MASTERMIND

The program contained in this directory allows you to play Mastermind* both in interactive mode and through an automatic simulation of gameplay run by the program.

The player, by running the interactive_player.py script, can specify the number of colors and the length of the sequence, which the program will generate randomly to then let the player guess it. The colors are represented by numbers - each number is a different color. After each player's attempt, the terminal receives a response from the judge, which gives the number of spots in the sequence that were hit and the number of colors the player gave that were out of place.

The auto_player.py script, meanwhile, accepts a sequence from the user for the computer to guess, using Donald Knuth's algorithm**.

## Table of Contents ##.

1. [Launching the program]
2. [Example 1 (interactive mode)]
3. [Example 2 (automatic mode)]
4. [Creators]

1. Launching the program

	!!! Python3.9 is required for proper operation of the program !!!

	1. in the terminal, change the directory to the one where this file is located by the command: cd directory path
	2. then giving the alias of the python interpreter, and after the space the name of the corresponding script, follow the prompts displayed in the terminal.
		(a) auto_player.py
			- provide a sequence (a sequence of numbers separated by spaces)
			- specify the number of colors considered (range of numbers) by typing the largest number from the range
		(b) interactive_player.py
			- specify the length of the sequence to be guessed
			- specify the number of colors considered (range of numbers) by typing the largest number from the range
			- try to guess the code by successive queries and analyzing the returned feedback
			- in order to terminate the program before guessing the hidden code, press {ctrl+D} (Linux) or {ctrl+Z then Enter} (Windows)

2 Example 1 (interactive mode)


	> cd /home/user/projects/Mastermind
	> python3 interactive_player.py

	Enter the number of colours:
	>6

	Enter the length of the sequence:
	>4

	Let's start the game!

	Enter your guess:
	>5 2 6 3
	You guessed 0 positions correctly
	among the missed ones, 2 of the colours you gave are on a wrong place
	Take another try :)

	Enter your guess:
	>6 4 2 4

	Nice! You have guessed the correct sequence in 2 tries!


3. Example 2 (interactive mode)

	> cd /home/user/projects/Mastermind
	> python3 auto_player.py

	Enter a sequence you want the computer to guess:
	>4 3 2 2

	Enter the number of colours considered:
	>6

	Computer is about to take up the challange!

	Round 1. Guess: (1, 1, 2, 2), Feedback: (2, 0), Number of possible answers left: 1296
	Round 2. Guess: (1, 2, 3, 4), Feedback: (0, 3), Number of possible answers left: 114
	Round 3. Guess: (1, 3, 2, 5), Feedback: (2, 0), Number of possible answers left: 16
	Round 4. Guess: (4, 3, 2, 2), Feedback: (4, 0), Number of possible answers left: 1
	Code found: (4, 3, 2, 2)
	Guesses: [(1, 1, 2, 2), (1, 2, 3, 4), (1, 3, 2, 5), (4, 3, 2, 2)]
	Done! Computer broke your code in 4 attempts.
	
4. Contributors

	Kacper Pietrzyk
	github:	https://github.com/AncientG7eek?tab=repositories
	email:	km.pietrzyk@stud.mimuw.edu.pl
