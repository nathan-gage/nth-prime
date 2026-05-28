CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -ffast-math -DNDEBUG -march=native -funroll-loops -flto -Wall -Wextra -pedantic
TARGET := nprime
SRC := src/nprime.cpp
WASM := nprime.wasm
WASM_SRC := web/nprime_wasm.rs
SITE_DIR := _site
WEB_ASSETS := index.html app.js styles.css

.PHONY: all clean test bench bench-scale profile-scale profile-hotpaths wasm pages

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

wasm: $(WASM)

$(WASM): $(WASM_SRC)
	rustc --target wasm32-unknown-unknown -O -C opt-level=3 -C lto=fat -C panic=abort --crate-type=cdylib $< -o $@

pages: $(WASM) $(WEB_ASSETS)
	rm -rf $(SITE_DIR)
	mkdir -p $(SITE_DIR)
	cp $(WEB_ASSETS) $(WASM) $(SITE_DIR)/
	touch $(SITE_DIR)/.nojekyll

test: $(TARGET)
	./$(TARGET) --self-test

bench: $(TARGET)
	./$(TARGET) --bench 1000000 10
	./$(TARGET) --bench 1000000000 10

bench-scale: $(TARGET)
	./$(TARGET) --bench 1000000000 10
	./$(TARGET) --bench 9999999999 3
	./$(TARGET) --bench 10000000000 3
	./$(TARGET) --bench 20000000000 1

profile-scale: $(TARGET)
	./$(TARGET) --profile 1000000
	./$(TARGET) --profile 1000000000
	./$(TARGET) --profile 9999999999

profile-hotpaths: $(TARGET)
	./$(TARGET) --profile-detail 1000000
	./$(TARGET) --profile-detail 200000000
	./$(TARGET) --profile-detail 1000000000
	./$(TARGET) --profile-detail 2000000000
	./$(TARGET) --profile-detail 5000000000
	./$(TARGET) --profile-detail 9999999999

clean:
	rm -f $(TARGET) $(WASM)
	rm -rf $(SITE_DIR)
