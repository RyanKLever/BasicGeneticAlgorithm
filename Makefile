macos:
	clang++ -std=c++17 -o run-main main.cpp
	clang++ -std=c++17 -o run-extra-features extra-features.cpp
gnu:
	g++ -std=c++17 -o run-main main.cpp
	g++ -std=c++17 -o run-extra-features extra-features.cpp
debug:
	g++ -std=c++17 -g -o run-main main.cpp
	g++ -std=c++17 -g -o run-extra-features extra-features.cpp
clean:
	rm -rf run-main*
	rm -rf run-extra-features*
