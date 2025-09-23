independents = {};
IC = 100;
expressions = 10000;

letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
letter_dict = {};
for (i = 0; i < Abs (letters); i+=1) {
    letter_dict + letters[i];
}

for (i = 0; i < IC; i+=1) {
    independents + Join ("", Random (letter_dict,1))[0][9];
    Eval (independents[i]);
    ^(independents[i]) = 1;
}

for (i = 0; i < expressions; i += 1) {
    
}




