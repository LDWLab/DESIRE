# Generated by Django 2.1.3 on 2019-10-10 16:23

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ribovision', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Mastertable',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('master_id', models.IntegerField()),
                ('active', models.IntegerField(blank=True, db_column='Active', null=True)),
                ('speciesname', models.CharField(blank=True, db_column='SpeciesName', max_length=24, null=True)),
                ('datasettype', models.CharField(blank=True, db_column='DataSetType', max_length=33, null=True)),
                ('structurename', models.CharField(blank=True, db_column='StructureName', max_length=14, null=True)),
                ('loadstring', models.CharField(blank=True, db_column='LoadString', max_length=29, null=True)),
                ('species_abr', models.CharField(blank=True, db_column='Species_Abr', max_length=5, null=True)),
            ],
            options={
                'db_table': 'MasterTable',
                'managed': False,
            },
        ),
    ]
