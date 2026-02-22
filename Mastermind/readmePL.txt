MASTERMIND

Zawarty w tym katalogu program umożliwia grę w Mastermind* zarówno w trybie interaktywnym, jak i poprzez automatyczną symulację rozgrywki prowadzoną przez program.

Gracz uruchamiając skrypt interactive_player.py można podać liczbę kolorów i długość sekwencji, którą program wygeneruje losowo, by następnie dać graczowi ją odgadnąć. Kolory reprezentowane są przez liczby - każda z nich to inny kolor. Po każdej próbie gracza, w terminalu otrzymuje on odpowiedź od sędziego, który podaje liczbę trafionych miejsc w sekwencji oraz liczbę kolorów które gracz podał nie na swoim miejscu.

Skrypt auto_player.py przyjmuje natomiast od użytkownika sekwencję, którą ma odgadnąć komputer, przy zastosowaniu algorytmu Donalda Knutha**.

*	https://pl.wikipedia.org/wiki/Mastermind_(gra_planszowa)
**	https://en.wikipedia.org/wiki/Mastermind_(board_game)#Worst_case:_Five-guess_algorithm


## Spis treści ##

1. [Uruchomienie programu]
2. [Przykład 1 (tryb interaktywny)]
3. [Przykład 2 (tryb automatyczny)]
4. [Twórcy]

1. Uruchomienie programu

	!!! Do poprawnej obsługi programu wymagany jest Python3.9 !!!

	1. W terminalu należy zmienić katalog na ten, w którym znajduje się ten plik poprzez komendę: cd ścieżka do katalogu
	2. Następnie podając alias interpretera pythona, a po spacji nazwę odpowiedniego skryptu, postępować zgodnie z wyświetlanymi w terminalu (w języku angielskim) komunikatami.
		a) auto_player.py
			- podać sekwencję (ciąg liczb oddzielonych spacjami)
			- podać liczbę rozpatrywanych kolorów (zakres liczb) wpisując największą liczbę z zakresu
		b) interactive_player.py
			- podać długość sekwencji do zgadnięcia
			- podać liczbę rozpatrywanych kolorów (zakres liczb) wpisując największą liczbę z zakresu
			- próbować odgadnąć kod poprzez kolejne zapytania i analizę zwracanego feedbacku
			- w celu zakończenia działania programu przed odgadnięciem ukrytego kodu, należy wcisnąć {ctrl+D} (Linux) lub {ctrl+Z następnie Enter} (Windows)

2. Przykład 1 (tryb interaktywny)

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


3. Przykład 2 (tryb automatyczny)

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
	
4. Twórcy

	Kacper Pietrzyk
	github:	https://github.com/AncientG7eek?tab=repositories
	email:	km.pietrzyk@stud.mimuw.edu.pl

	