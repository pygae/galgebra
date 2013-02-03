# Configuration file for ipython-qtconsole.

c = get_config()

#------------------------------------------------------------------------------
# IPythonQtConsoleApp configuration
#------------------------------------------------------------------------------

# IPythonQtConsoleApp will inherit config from: BaseIPythonApplication,
# Application, IPythonConsoleApp

# The IPython profile to use.
# c.IPythonQtConsoleApp.profile = u'default'

# set the stdin (XREQ) port [default: random]
# c.IPythonQtConsoleApp.stdin_port = 0

# Set the log level by value or name.
# c.IPythonQtConsoleApp.log_level = 30

# Path to the ssh key to use for logging in to the ssh server.
# c.IPythonQtConsoleApp.sshkey = ''

# Use a plaintext widget instead of rich text (plain can't print/save).
# c.IPythonQtConsoleApp.plain = False

# Set the kernel's IP address [default localhost]. If the IP address is
# something other than localhost, then Consoles on other machines will be able
# to connect to the Kernel, so be careful!
# c.IPythonQtConsoleApp.ip = '127.0.0.1'

# JSON file in which to store connection info [default: kernel-<pid>.json]
# 
# This file will contain the IP, ports, and authentication key needed to connect
# clients to this kernel. By default, this file will be created in the security-
# dir of the current profile, but can be specified by absolute path.
# c.IPythonQtConsoleApp.connection_file = ''

# Connect to an already running kernel
# c.IPythonQtConsoleApp.existing = ''

# Create a massive crash report when IPython enconters what may be an internal
# error.  The default is to append a short message to the usual traceback
# c.IPythonQtConsoleApp.verbose_crash = False

# path to a custom CSS stylesheet
# c.IPythonQtConsoleApp.stylesheet = ''

# set the heartbeat port [default: random]
# c.IPythonQtConsoleApp.hb_port = 0

# The name of the IPython directory. This directory is used for logging
# configuration (through profiles), history storage, etc. The default is usually
# $HOME/.ipython. This options can also be specified through the environment
# variable IPYTHON_DIR.
# c.IPythonQtConsoleApp.ipython_dir = u'/home/brombo/.config/ipython'

# The SSH server to use to connect to the kernel.
# c.IPythonQtConsoleApp.sshserver = ''

# Use a pure Python kernel instead of an IPython kernel.
# c.IPythonQtConsoleApp.pure = False

# Set to display confirmation dialog on exit. You can always use 'exit' or
# 'quit', to force a direct exit without any confirmation.
# c.IPythonQtConsoleApp.confirm_exit = True

# Whether to install the default config files into the profile dir. If a new
# profile is being created, and IPython contains config files for that profile,
# then they will be staged into the new directory.  Otherwise, default config
# files will be automatically generated.
# c.IPythonQtConsoleApp.copy_config_files = False

# set the shell (XREP) port [default: random]
# c.IPythonQtConsoleApp.shell_port = 0

# Whether to create profile dir if it doesn't exist
# c.IPythonQtConsoleApp.auto_create = False

# Whether to overwrite existing config files when copying
# c.IPythonQtConsoleApp.overwrite = False

# set the iopub (PUB) port [default: random]
# c.IPythonQtConsoleApp.iopub_port = 0

#------------------------------------------------------------------------------
# IPKernelApp configuration
#------------------------------------------------------------------------------

# IPython: an enhanced interactive Python shell.

# IPKernelApp will inherit config from: KernelApp, BaseIPythonApplication,
# Application, InteractiveShellApp

# The importstring for the DisplayHook factory
# c.IPKernelApp.displayhook_class = 'IPython.zmq.displayhook.ZMQDisplayHook'

# Set the IP or interface on which the kernel will listen.
# c.IPKernelApp.ip = '127.0.0.1'

# Pre-load matplotlib and numpy for interactive use, selecting a particular
# matplotlib backend and loop integration.
# c.IPKernelApp.pylab = None

# Create a massive crash report when IPython enconters what may be an internal
# error.  The default is to append a short message to the usual traceback
# c.IPKernelApp.verbose_crash = False

# set the shell (XREP) port [default: random]
# c.IPKernelApp.shell_port = 0

# Whether to overwrite existing config files when copying
# c.IPKernelApp.overwrite = False

# Execute the given command string.
# c.IPKernelApp.code_to_run = ''

# set the stdin (XREQ) port [default: random]
# c.IPKernelApp.stdin_port = 0

# Set the log level by value or name.
# c.IPKernelApp.log_level = 30

# lines of code to run at IPython startup.
# c.IPKernelApp.exec_lines = []

# The importstring for the OutStream factory
# c.IPKernelApp.outstream_class = 'IPython.zmq.iostream.OutStream'

# Whether to create profile dir if it doesn't exist
# c.IPKernelApp.auto_create = False

# set the heartbeat port [default: random]
# c.IPKernelApp.hb_port = 0

# redirect stdout to the null device
# c.IPKernelApp.no_stdout = False

# dotted module name of an IPython extension to load.
# c.IPKernelApp.extra_extension = ''

# A file to be run
# c.IPKernelApp.file_to_run = ''

# The IPython profile to use.
# c.IPKernelApp.profile = u'default'

# 
# c.IPKernelApp.parent_appname = u''

# kill this process if its parent dies.  On Windows, the argument specifies the
# HANDLE of the parent process, otherwise it is simply boolean.
# c.IPKernelApp.parent = 0

# JSON file in which to store connection info [default: kernel-<pid>.json]
# 
# This file will contain the IP, ports, and authentication key needed to connect
# clients to this kernel. By default, this file will be created in the security-
# dir of the current profile, but can be specified by absolute path.
# c.IPKernelApp.connection_file = ''

# If true, an 'import *' is done from numpy and pylab, when using pylab
# c.IPKernelApp.pylab_import_all = True

# The name of the IPython directory. This directory is used for logging
# configuration (through profiles), history storage, etc. The default is usually
# $HOME/.ipython. This options can also be specified through the environment
# variable IPYTHON_DIR.
# c.IPKernelApp.ipython_dir = u'/home/brombo/.config/ipython'

# ONLY USED ON WINDOWS Interrupt this process when the parent is signalled.
# c.IPKernelApp.interrupt = 0

# Whether to install the default config files into the profile dir. If a new
# profile is being created, and IPython contains config files for that profile,
# then they will be staged into the new directory.  Otherwise, default config
# files will be automatically generated.
# c.IPKernelApp.copy_config_files = False

# List of files to run at IPython startup.
# c.IPKernelApp.exec_files = []

# A list of dotted module names of IPython extensions to load.
# c.IPKernelApp.extensions = []

# redirect stderr to the null device
# c.IPKernelApp.no_stderr = False

# set the iopub (PUB) port [default: random]
# c.IPKernelApp.iopub_port = 0

#------------------------------------------------------------------------------
# IPythonWidget configuration
#------------------------------------------------------------------------------

# A FrontendWidget for an IPython kernel.

# IPythonWidget will inherit config from: FrontendWidget, HistoryConsoleWidget,
# ConsoleWidget

# 
# c.IPythonWidget.input_sep = '\n'

# The type of underlying text widget to use. Valid values are 'plain', which
# specifies a QPlainTextEdit, and 'rich', which specifies a QTextEdit.
# c.IPythonWidget.kind = 'plain'

# 
# c.IPythonWidget.output_sep2 = ''

# The font size. If unconfigured, Qt will be entrusted with the size of the
# font.
# c.IPythonWidget.font_size = 0

# Whether to draw information calltips on open-parentheses.
# c.IPythonWidget.enable_calltips = True

# Use a list widget instead of plain text output for tab completion.
# c.IPythonWidget.gui_completion = False

# 
# c.IPythonWidget.in_prompt = 'In [<span class="in-prompt-number">%i</span>]: '

# Whether to process ANSI escape codes.
# c.IPythonWidget.ansi_codes = True

# The editor command to use when a specific line number is requested. The string
# should contain two format specifiers: {line} and {filename}. If this parameter
# is not specified, the line number option to the %edit magic will be ignored.
# c.IPythonWidget.editor_line = u''

# A CSS stylesheet. The stylesheet can contain classes for:
#     1. Qt: QPlainTextEdit, QFrame, QWidget, etc
#     2. Pygments: .c, .k, .o, etc. (see PygmentsHighlighter)
#     3. IPython: .error, .in-prompt, .out-prompt, etc
# c.IPythonWidget.style_sheet = u''

# The type of paging to use. Valid values are:
# 
#     'inside' : The widget pages like a traditional terminal.
#     'hsplit' : When paging is requested, the widget is split
#                horizontally. The top pane contains the console, and the
#                bottom pane contains the paged text.
#     'vsplit' : Similar to 'hsplit', except that a vertical splitter
#                used.
#     'custom' : No action is taken by the widget beyond emitting a
#                'custom_page_requested(str)' signal.
#     'none'   : The text is written directly to the console.
# c.IPythonWidget.paging = 'inside'

# 
# c.IPythonWidget.output_sep = ''

# A command for invoking a system text editor. If the string contains a
# {filename} format specifier, it will be used. Otherwise, the filename will be
# appended to the end the command.
# c.IPythonWidget.editor = ''

# If not empty, use this Pygments style for syntax highlighting. Otherwise, the
# style sheet is queried for Pygments style information.
# c.IPythonWidget.syntax_style = u''

# The maximum number of lines of text before truncation. Specifying a non-
# positive number disables text truncation (not recommended).
# c.IPythonWidget.buffer_size = 500

# 
# c.IPythonWidget.history_lock = False

# The font family to use for the console. On OSX this defaults to Monaco, on
# Windows the default is Consolas with fallback of Courier, and on other
# platforms the default is Monospace.
# c.IPythonWidget.font_family = u''

# 
# c.IPythonWidget.out_prompt = 'Out[<span class="out-prompt-number">%i</span>]: '

#------------------------------------------------------------------------------
# ZMQInteractiveShell configuration
#------------------------------------------------------------------------------

# A subclass of InteractiveShell for ZMQ.

# ZMQInteractiveShell will inherit config from: InteractiveShell

# Use colors for displaying information about objects. Because this information
# is passed through a pager (like 'less'), and some pagers get confused with
# color codes, this capability can be turned off.
# c.ZMQInteractiveShell.color_info = True

# 
# c.ZMQInteractiveShell.history_length = 10000

# Don't call post-execute functions that have failed in the past.
# c.ZMQInteractiveShell.disable_failing_post_execute = False

# Show rewritten input, e.g. for autocall.
# c.ZMQInteractiveShell.show_rewritten_input = True

# Set the color scheme (NoColor, Linux, or LightBG).
# c.ZMQInteractiveShell.colors = 'Linux'

# 
# c.ZMQInteractiveShell.separate_in = '\n'

# Deprecated, use PromptManager.in2_template
# c.ZMQInteractiveShell.prompt_in2 = '   .\\D.: '

# 
# c.ZMQInteractiveShell.separate_out = ''

# Deprecated, use PromptManager.in_template
# c.ZMQInteractiveShell.prompt_in1 = 'In [\\#]: '

# Enable deep (recursive) reloading by default. IPython can use the deep_reload
# module which reloads changes in modules recursively (it replaces the reload()
# function, so you don't need to change anything to use it). deep_reload()
# forces a full reload of modules whose code may have changed, which the default
# reload() function does not.  When deep_reload is off, IPython will use the
# normal reload(), but deep_reload will still be available as dreload().
# c.ZMQInteractiveShell.deep_reload = False

# Make IPython automatically call any callable object even if you didn't type
# explicit parentheses. For example, 'str 43' becomes 'str(43)' automatically.
# The value can be '0' to disable the feature, '1' for 'smart' autocall, where
# it is not applied if there are no more arguments on the line, and '2' for
# 'full' autocall, where all callable objects are automatically called (even if
# no arguments are present).
# c.ZMQInteractiveShell.autocall = 0

# 
# c.ZMQInteractiveShell.separate_out2 = ''

# Deprecated, use PromptManager.justify
# c.ZMQInteractiveShell.prompts_pad_left = True

# 
# c.ZMQInteractiveShell.readline_parse_and_bind = ['tab: complete', '"\\C-l": clear-screen', 'set show-all-if-ambiguous on', '"\\C-o": tab-insert', '"\\C-r": reverse-search-history', '"\\C-s": forward-search-history', '"\\C-p": history-search-backward', '"\\C-n": history-search-forward', '"\\e[A": history-search-backward', '"\\e[B": history-search-forward', '"\\C-k": kill-line', '"\\C-u": unix-line-discard']

# Enable magic commands to be called without the leading %.
# c.ZMQInteractiveShell.automagic = True

# 
# c.ZMQInteractiveShell.debug = False

# 
# c.ZMQInteractiveShell.object_info_string_level = 0

# 
# c.ZMQInteractiveShell.ipython_dir = ''

# 
# c.ZMQInteractiveShell.readline_remove_delims = '-/~'

# Start logging to the default log file.
# c.ZMQInteractiveShell.logstart = False

# The name of the logfile to use.
# c.ZMQInteractiveShell.logfile = ''

# 
# c.ZMQInteractiveShell.wildcards_case_sensitive = True

# Save multi-line entries as one entry in readline history
# c.ZMQInteractiveShell.multiline_history = True

# Start logging to the given file in append mode.
# c.ZMQInteractiveShell.logappend = ''

# 
# c.ZMQInteractiveShell.xmode = 'Context'

# 
# c.ZMQInteractiveShell.quiet = False

# Deprecated, use PromptManager.out_template
# c.ZMQInteractiveShell.prompt_out = 'Out[\\#]: '

# Set the size of the output cache.  The default is 1000, you can change it
# permanently in your config file.  Setting it to 0 completely disables the
# caching system, and the minimum value accepted is 20 (if you provide a value
# less than 20, it is reset to 0 and a warning is issued).  This limit is
# defined because otherwise you'll spend more time re-flushing a too small cache
# than working
# c.ZMQInteractiveShell.cache_size = 1000

# Automatically call the pdb debugger after every exception.
# c.ZMQInteractiveShell.pdb = False

#------------------------------------------------------------------------------
# ProfileDir configuration
#------------------------------------------------------------------------------

# An object to manage the profile directory and its resources.
# 
# The profile directory is used by all IPython applications, to manage
# configuration, logging and security.
# 
# This object knows how to find, create and manage these directories. This
# should be used by any code that wants to handle profiles.

# Set the profile location directly. This overrides the logic used by the
# `profile` option.
# c.ProfileDir.location = u''

#------------------------------------------------------------------------------
# Session configuration
#------------------------------------------------------------------------------

# Object for handling serialization and sending of messages.
# 
# The Session object handles building messages and sending them with ZMQ sockets
# or ZMQStream objects.  Objects can communicate with each other over the
# network via Session objects, and only need to work with the dict-based IPython
# message spec. The Session will handle serialization/deserialization, security,
# and metadata.
# 
# Sessions support configurable serialiization via packer/unpacker traits, and
# signing with HMAC digests via the key/keyfile traits.
# 
# Parameters ----------
# 
# debug : bool
#     whether to trigger extra debugging statements
# packer/unpacker : str : 'json', 'pickle' or import_string
#     importstrings for methods to serialize message parts.  If just
#     'json' or 'pickle', predefined JSON and pickle packers will be used.
#     Otherwise, the entire importstring must be used.
# 
#     The functions must accept at least valid JSON input, and output *bytes*.
# 
#     For example, to use msgpack:
#     packer = 'msgpack.packb', unpacker='msgpack.unpackb'
# pack/unpack : callables
#     You can also set the pack/unpack callables for serialization directly.
# session : bytes
#     the ID of this Session object.  The default is to generate a new UUID.
# username : unicode
#     username added to message headers.  The default is to ask the OS.
# key : bytes
#     The key used to initialize an HMAC signature.  If unset, messages
#     will not be signed or checked.
# keyfile : filepath
#     The file containing a key.  If this is set, `key` will be initialized
#     to the contents of the file.

# Username for the Session. Default is your system username.
# c.Session.username = 'brombo'

# The name of the packer for serializing messages. Should be one of 'json',
# 'pickle', or an import name for a custom callable serializer.
# c.Session.packer = 'json'

# The UUID identifying this session.
# c.Session.session = u''

# execution key, for extra authentication.
# c.Session.key = ''

# Debug output in the Session
# c.Session.debug = False

# The name of the unpacker for unserializing messages. Only used with custom
# functions for `packer`.
# c.Session.unpacker = 'json'

# path to file containing execution key.
# c.Session.keyfile = ''
