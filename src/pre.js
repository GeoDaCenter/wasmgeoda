const base64 = require('base64-js');
const base64EncodedBinary = require('./jsgeoda-wasm');

// The emscripten's Module object.
// See http://kripken.github.io/emscripten-site/docs/api_reference/module.html for details.
var Module = {};

Module.noInitialRun = true;
Module.wasmBinary = base64.toByteArray(base64EncodedBinary);
