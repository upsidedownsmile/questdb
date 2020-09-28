```
const net = require('net');
const { Crypto } = require("node-webcrypto-ossl");
const crypto = new Crypto();
const base64 = require('base64-arraybuffer');

let jwk = {
  "kid": "testUser1",
  "kty": "EC",
  "d": "5UjEMuA0Pj5pjK8a-fa24dyIf-Es5mYny3oE_Wmus48",
  "crv": "P-256",
  "x": "fLKYEaoEb9lrn3nkwLDA-M_xnuFOdSt9y0Z7_vWSHLU",
  "y": "Dt5tbS1dEDMSYfym3fgMv0B99szno-dFc1rYF9t0aac"
};

async function sendTestMeasurement() {
 let apiKey = await crypto.subtle.importKey("jwk", jwk, {name: "ECDSA", namedCurve: "P-256"}, true, ["sign", "verify"]);

 let client = new net.Socket();
 let challenge;
 client.on('data', async function(data) {
  if (typeof challenge === 'undefined') {
   challenge = data;
  } else {
   challenge = Buffer.concat([challenge, data]);
  }

  // Check for trailing \n which ends the challenge
  if (challenge.slice(challenge.length-1).readInt8() == 10) {
   let rawSignature = await crypto.subtle.sign({name: "ECDSA", hash: "SHA-256"}, apiKey, challenge.slice(0, challenge.length-1));
   let signature = base64.encode(rawSignature);
   client.write(signature+'\n');

   // Send measurement
   let measurement = 'test,location=us-midwest temperature=82.4 '+(new Date() * 1000000);
   console.log('measurement: ' + measurement);
   client.write(measurement+'\n');
   client.destroy();
  }
 });

 client.on('close', function() {
  console.log('Connection closed');
 });

 client.connect(9009, 'ilp.savvy.questdb.net', function() {
  client.write(jwk.kid+'\n');
 });
}

sendTestMeasurement();
```