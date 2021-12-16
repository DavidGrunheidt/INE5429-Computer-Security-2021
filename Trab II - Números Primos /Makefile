.PHONY: start

start:
	gcc main.c -lgmp -I /opt/homebrew/Cellar/gmp/6.2.1_1/include -L /opt/homebrew/Cellar/gmp/6.2.1_1/lib -o main
	chmod +X main
	./main | tee output

clean:
	git clean -fdx