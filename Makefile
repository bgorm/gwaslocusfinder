.DEFAULT_GOAL := all
CXX ?= c++
CXXFLAGS ?= -O3 -g
WITH_PGEN ?= 1
BUILD_DIR ?= build
HTSLIB ?= vendor/htslib-1.25
PLUGIN_HTSLIB ?= $(HTSLIB)
PGENLIB ?= vendor/plink-ng-4544f92e374e784a6ffd97fa88239962b81f2b62/2.0/include
BCFTOOLS ?= bcftools
BCFTOOLS_VERSION ?= $(shell $(BCFTOOLS) --version-only 2>/dev/null | cut -d+ -f1)
PLUGIN_HTSLIB_VERSION ?= $(shell cd $(PLUGIN_HTSLIB) && sh version.sh 2>/dev/null)
PYTHON ?= python3
PKG_CONFIG ?= pkg-config
PLATFORM ?= $(shell uname -s)
PREFIX ?= /usr/local
plugindir ?= $(PREFIX)/libexec/bcftools

ifneq ($(WITH_PGEN),0)
ifneq ($(WITH_PGEN),1)
$(error WITH_PGEN must be 0 or 1)
endif
endif

# pkg-config helps with non-system prefixes (for example Homebrew). All flags
# can be overridden; CPPFLAGS, CXXFLAGS, LDFLAGS and LDLIBS are also respected.
COMPRESSION_CPPFLAGS ?= $(shell $(PKG_CONFIG) --cflags zlib libzstd 2>/dev/null)
COMPRESSION_LIBS ?= $(shell $(PKG_CONFIG) --libs zlib libzstd 2>/dev/null || echo -lz -lzstd)
PROJECT_CPPFLAGS = -Isrc -DLF_WITH_PGEN=$(WITH_PGEN) $(COMPRESSION_CPPFLAGS)
HEADERS = src/core.hpp src/io.hpp src/storage.hpp src/reference.hpp src/entry.hpp
OBJDIR = $(BUILD_DIR)/plugin-pgen$(WITH_PGEN)
PG_NAMES =
ifeq ($(WITH_PGEN),1)
PROJECT_CPPFLAGS += -I$(PGENLIB)
PG_NAMES = plink2_base plink2_bits plink2_simd pgenlib_misc pgenlib_read
endif
PG_PIC_OBJECTS = $(addprefix $(OBJDIR)/,$(addsuffix .o,$(PG_NAMES)))
PLUGIN_OBJECTS = $(OBJDIR)/engine.o $(OBJDIR)/wrapper.o $(PG_PIC_OBJECTS)
PLUGIN = $(BUILD_DIR)/gwaslocusfinder.so

ifeq ($(PLATFORM),Darwin)
BCFTOOLS_PATH := $(shell command -v $(BCFTOOLS))
PLUGIN_LDFLAGS ?= -bundle -bundle_loader "$(BCFTOOLS_PATH)" -Wl,-undefined,dynamic_lookup
else
PLUGIN_LDFLAGS ?= -shared
endif

.PHONY: all plugin test test-storage install export clean FORCE
all: plugin
plugin: $(PLUGIN)
$(OBJDIR):
	mkdir -p "$@"
$(OBJDIR)/%.o: $(PGENLIB)/%.cc | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CXXFLAGS) -std=c++17 -DNDEBUG -fPIC -fvisibility=hidden -c $< -o $@
$(OBJDIR)/engine.o: src/engine.cpp $(HEADERS) $(PLUGIN_HTSLIB)/htslib/hts.h | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) -I$(PLUGIN_HTSLIB) $(CXXFLAGS) -std=c++17 -Wall -Wextra -Wpedantic -fPIC -fvisibility=hidden -c $< -o $@
$(OBJDIR)/wrapper.o: src/plugin.cpp src/entry.hpp $(PLUGIN_HTSLIB)/htslib/hts.h $(wildcard $(PLUGIN_HTSLIB)/version.sh) | $(OBJDIR)
	@test -n "$(BCFTOOLS_VERSION)" || { echo "bcftools is required; set BCFTOOLS=/path/to/bcftools" >&2; exit 1; }
	@test -n "$(PLUGIN_HTSLIB_VERSION)" || { echo "Set PLUGIN_HTSLIB_VERSION when using installed HTSlib headers" >&2; exit 1; }
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) -I$(PLUGIN_HTSLIB) $(CXXFLAGS) -std=c++17 -fPIC -fvisibility=hidden \
	  -DLF_BCFTOOLS_VERSION='"$(BCFTOOLS_VERSION)"' -DLF_HTSLIB_VERSION='"$(PLUGIN_HTSLIB_VERSION)"' -c $< -o $@
# Force a relink when selecting an already-built configuration: an older set
# of PGEN-enabled objects must not leave a newer VCF-only plugin in place.
# HTSlib and bcftools symbols resolve from the host; no second HTSlib is linked.
$(PLUGIN): $(PLUGIN_OBJECTS) FORCE
	$(CXX) $(CXXFLAGS) $(PLUGIN_LDFLAGS) $(LDFLAGS) $(PLUGIN_OBJECTS) $(COMPRESSION_LIBS) -lm -lpthread $(LDLIBS) -o $@
$(BUILD_DIR)/test-storage: tests/test_storage.cpp src/storage.hpp src/io.hpp src/core.hpp | $(OBJDIR)
	$(CXX) $(CPPFLAGS) $(PROJECT_CPPFLAGS) -I$(PLUGIN_HTSLIB) $(CXXFLAGS) -std=c++17 -Wall -Wextra -Wpedantic $< -o $@
test-storage: $(BUILD_DIR)/test-storage
	$(BUILD_DIR)/test-storage
test: plugin test-storage
	LF_PLUGIN="$(abspath $(PLUGIN))" LF_EXPECT_PGEN="$(WITH_PGEN)" BCFTOOLS="$(BCFTOOLS)" $(PYTHON) -B -m unittest discover -s tests -p test_plugin.py -v
install: plugin
	mkdir -p "$(DESTDIR)$(plugindir)"
	install -m 644 $(PLUGIN) "$(DESTDIR)$(plugindir)/gwaslocusfinder.so"
export:
	@test -n "$(EXPORT_DIR)" || { echo "Usage: make export EXPORT_DIR=/path/to/new-source-directory" >&2; exit 1; }
	bash scripts/export-plugin.sh "$(EXPORT_DIR)"
clean:
	rm -f $(BUILD_DIR)/test-storage $(BUILD_DIR)/gwaslocusfinder.so $(BUILD_DIR)/plugin-pgen0/*.o $(BUILD_DIR)/plugin-pgen1/*.o
FORCE:

# Local development targets are intentionally omitted from public exports.
-include Makefile.internal
