from setuptools import setup

setup(
    name='HCI-FITS',
    version='0.1',

    description='Module dedicated to the creation of high-level science products for high-contrast imaging observations (HCI-FITS format)',
    url='https://github.com/avigan/HCI-FITS',
    author='Arthur Vigan',
    author_email='arthur.vigan@lam.fr',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Professional Astronomers',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    keywords='image processing, data analysis, high-contrast imaging',
    packages=['hcifits'],
    install_requires=[
        'numpy', 'scipy', 'astropy'
    ],
    zip_safe=False
)
