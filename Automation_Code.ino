
/*
  Cleaner Pattern State-Machine
  Reads digital inputs on pins 0,1,2,5,9 and analog on A1,A2.
  Drives outputs on pins 3,4,8 based on a simple 4-state machine.
*/

int I = 0;  // State variable

void setup() {
  // Outputs
  pinMode(3, OUTPUT);
  pinMode(4, OUTPUT);
  pinMode(8, OUTPUT);
  // Inputs
  pinMode(0, INPUT);
  pinMode(1, INPUT);
  pinMode(2, INPUT);
  pinMode(5, INPUT);
  pinMode(9, INPUT);
}

void loop() {
  // Re-read inputs each pass
  bool a = digitalRead(0);
  bool b = digitalRead(1);
  bool c = digitalRead(2);
  bool e = digitalRead(5);
  bool i = digitalRead(9);
  int f = analogRead(A1);
  int g = analogRead(A2);

  // State-machine
  switch (I) {
    case 0:  // Idle
      if (a) {
        I = 1;
        digitalWrite(3, LOW);
        digitalWrite(4, HIGH);
      }
      else if (b) {
        I = 2;
        digitalWrite(3, HIGH);
        digitalWrite(4, HIGH);
      }
      break;

    case 1:  // Active after 'a'
    case 2:  // Active after 'b'
      if (c || e || i) {
        I = 3;
        digitalWrite(4, LOW);
      }
      break;

    case 3:  // Done
      // No further action
      break;
  }

  // Analog-triggered outputs
  if (I == 0 && f > 10) {
    digitalWrite(8, HIGH);
  }
  if ((I == 1 || I == 2) && g > 15) {
    digitalWrite(4, LOW);
  }
}
