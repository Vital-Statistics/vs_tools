#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:44:08 2024

@author: rudy
"""

# Function to encrypt data
def encrypt_file(file_path, key=None):
    import os
    import base64
    from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
    from cryptography.hazmat.primitives import padding
    from cryptography.hazmat.backends import default_backend
    
    if key is None:
        key=base64.b64decode(os.environ['ENCRYPTION_KEY'])
    with open(file_path, 'rb') as file:
        file_data = file.read()
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())
    encryptor = cipher.encryptor()
    padder = padding.PKCS7(algorithms.AES.block_size).padder()
    padded_data = padder.update(file_data) + padder.finalize()
    encrypted_data = encryptor.update(padded_data) + encryptor.finalize()
    with open(file_path + '.enc', 'wb') as encrypted_file:
        encrypted_file.write(iv + encrypted_data)

